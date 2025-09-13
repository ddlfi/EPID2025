#include "witness_prove.hpp"
#include "epid.h"
#include "aes_prove.h"
#include "parameters.h"
#include "utils.h"
#include "compat.h"

// SRL prove function using signature_t
static void srl_prove_128(const uint8_t* w, const bf128_t* bf_v,
                          bf128_t* a_0_vec, bf128_t* a_1_vec, bf128_t* a_2_vec,
                          unsigned int lambda, uint8_t* k, bf128_t* vk,
                          const std::vector<signature_t>& srl_list) {
    unsigned int offset = 0;
    unsigned int a_vec_offset = 0;

    for (auto& sig : srl_list) {
        // AES prove: f(sk_i, r_j) with in_flag=0, out_flag=1
        // r_j is public input, aes_out is secret
        aes_prove_128(w + offset / 8, bf_v + offset, sig.r.data(),
                      w + offset / 8 + 9 * 128 / 8, a_0_vec + a_vec_offset,
                      a_1_vec + a_vec_offset, a_2_vec + a_vec_offset, 0, k, vk, 0, 1);

        offset += AES_128_2;  // AES witness size (in bits)
        a_vec_offset += 320;

        // OR prove: prove cumulative OR of d bits
        or_prove(w + (offset - 128) / 8, bf_v + offset - 128,
                 a_0_vec + a_vec_offset, a_1_vec + a_vec_offset, a_2_vec + a_vec_offset, lambda,
                 sig.r.data(), sig.t.data());

        offset += 128;  // OR witness size (in bits)
        a_vec_offset += (lambda - 1);
    }
}

// SRL verify function using signature_t
static void srl_verify_128(const bf128_t* bf_q, const uint8_t* delta,
                           bf128_t* b_vec, unsigned int lambda, bf128_t* qk,
                           const std::vector<signature_t>& srl_list) {
    unsigned int witness_offset = 0;
    unsigned int b_vec_offset = 0;

    for (auto& sig : srl_list) {
        // AES verify: f(sk_i, r_j) with in_flag=0, out_flag=1
        // r_j is public input, aes_out is secret
        aes_verify_128(bf_q + witness_offset, delta, sig.r.data(), nullptr,
                       b_vec + b_vec_offset, false, qk, false, true);

        witness_offset += AES_128_2;
        b_vec_offset += 320;

        // OR verify: verify cumulative OR of d bits
        or_verify(bf_q + witness_offset - 128, delta, b_vec + b_vec_offset,
                  lambda, sig.r.data(), sig.t.data());

        witness_offset += 128;         // OR witness size (in bf128_t units)
        b_vec_offset += (lambda - 1);  // OR constraints count
    }
}


void witness_prove(const uint8_t* w, const bf128_t* bf_v, bf128_t* a_0_vec,
                   bf128_t* a_1_vec, bf128_t* a_2_vec, uint8_t* chall_2, unsigned int lambda,
                   uint8_t* t, uint8_t* r, unsigned int height,
                   unsigned int ell_hat, const std::vector<signature_t>& srl) {
    unsigned int index = 0;
    uint8_t* k = (uint8_t*)malloc((10 + 1) * 128 / 8);
    bf128_t* vk = (bf128_t*)faest_aligned_alloc(
        BF128_ALIGN, sizeof(bf128_t) * ((10 + 1) * 128));

    // t = f(sk_i, r)
    aes_prove_128(w, bf_v, r, t, a_0_vec, a_1_vec, a_2_vec, true, k, vk, false, false);

    index += 128 + 10 * 32 + 9 * 128;
    // t_i_join = f(sk_i,c_i)
    aes_prove_128(w + index / 8, bf_v + index, w + (index + 9 * 128 + 128) / 8,
                  w + (index + 9 * 128) / 8, a_0_vec + 400, a_1_vec + 400, a_2_vec + 400,
                  false, k, vk, true, true);
    index += 9 * 128 + 128;  // states + aes_out

    // x_i = f(c_i,t_i_join)
    leave_prove_128(w + (index - 128) / 8, bf_v + (index - 128), a_0_vec + 720,
                    a_1_vec + 720, a_2_vec + 720);

    index += 128 + 10 * 32 + 9 * 128;

    // Merkle tree prove
    merkle_tree_prove_128(w + index / 8, bf_v + index, a_0_vec + 1120,
                          a_1_vec + 1120, a_2_vec + 1120, height);

    index += (128 * 3 + 128 + 10 * 32 + 128 * 9) * height + 128;

    // SRL prove
    if (!srl.empty()) {
        srl_prove_128(w + index / 8, bf_v + index, a_0_vec + 1120 + 400 * height,
                      a_1_vec + 1120 + 400 * height, a_2_vec + 1120 + 400 * height, lambda, k, vk, srl);
    }

    faest_aligned_free(vk);
    free(k);
}

void witness_verify(const bf128_t* bf_q, const uint8_t* delta, bf128_t* b_vec,
                    unsigned int lambda, const uint8_t* t, const uint8_t* r,
                    unsigned int height, unsigned int ell_hat,
                    const std::vector<signature_t>& srl) {
    unsigned int index = 0;
    bf128_t* qk = (bf128_t*)faest_aligned_alloc(
        BF128_ALIGN, sizeof(bf128_t) * ((10 + 1) * 128));

    aes_verify_128(bf_q, delta, r, t, b_vec, true, qk, false, false);
    index += 128 + 10 * 32 + 9 * 128;

    aes_verify_128(bf_q + index, delta, nullptr, nullptr, b_vec + 400, false,
                   qk, true, true);
    index += 9 * 128 + 128;

    leave_verify_128(bf_q + (index - 128), delta, b_vec + 720);

    index += 128 + 10 * 32 + 9 * 128;

    merkle_tree_verify_128(bf_q + index, delta, b_vec + 1120, height);

    index += (128 * 3 + 128 + 10 * 32 + 128 * 9) * height + 128;

    // SRL verify
    if (!srl.empty()) {
        srl_verify_128(bf_q + index, delta, b_vec + 1120 + 400 * height, lambda,
                       qk, srl);
    }

    faest_aligned_free(qk);
}


void or_prove(const uint8_t* w, const bf128_t* bf_v, bf128_t* a_0_vec,
              bf128_t* a_1_vec, bf128_t* a_2_vec, unsigned int lambda, const uint8_t* r_j,
              const uint8_t* t_j) {
    unsigned int n = lambda;

    uint8_t temp_d[lambda / 8];
    xor_u8_array(w, r_j, temp_d, lambda / 8);
    xor_u8_array(temp_d, t_j, temp_d, lambda / 8);

    uint8_t u1 = ptr_get_bit(temp_d, 0);
    uint8_t u2 = ptr_get_bit(temp_d, 1);
    bf128_t v1 = bf_v[0];
    bf128_t v2 = bf_v[1];
    bf128_t v_or = bf_v[n];
    bf128_t v3 = bf128_add(v_or, bf128_add(v1, v2));

    a_0_vec[0] = bf128_mul(v1, v2);
    a_1_vec[0] =
        bf128_add(bf128_add(bf128_mul_bit(v2, u1), bf128_mul_bit(v1, u2)), v3);

    for (unsigned int i = 1; i < n - 1; i++) {
        uint8_t u1 = ptr_get_bit(w, n + i - 1);
        uint8_t u2 = ptr_get_bit(temp_d, i + 1);

        bf128_t v1 = bf_v[n + i - 1];
        bf128_t v2 = bf_v[i + 1];
        bf128_t v_or = bf_v[n + i];
        bf128_t v3 = bf128_add(v_or, bf128_add(v1, v2));

        a_0_vec[i] = bf128_mul(v1, v2);
        a_1_vec[i] = bf128_add(
            bf128_add(bf128_mul_bit(v2, u1), bf128_mul_bit(v1, u2)), v3);
    }
}

void or_verify(const bf128_t* bf_q, const uint8_t* delta, bf128_t* b_vec,
               unsigned int lambda, const uint8_t* r_j, const uint8_t* t_j) {
    unsigned int n = lambda;
    const bf128_t bf_delta = bf128_load(delta);

    bf128_t* temp_q = (bf128_t*)malloc(128 * sizeof(bf128_t));
    uint8_t rt_xor[lambda / 8];
    xor_u8_array(r_j, t_j, rt_xor, lambda / 8);

    for (unsigned int bit_idx = 0; bit_idx < lambda; bit_idx++) {
        temp_q[bit_idx] = bf_q[bit_idx];
        if (bit_idx < lambda) {
            uint8_t rt_bit = ptr_get_bit(rt_xor, bit_idx);
            bf128_t adjustment = bf128_mul_bit(bf_delta, rt_bit);
            temp_q[bit_idx] = bf128_add(bf_q[bit_idx], adjustment);
        }
    }

    bf128_t q1 = temp_q[0];
    bf128_t q2 = temp_q[1];
    bf128_t q_or = bf_q[n];
    bf128_t q3 = bf128_add(q_or, bf128_add(q1, q2));

    b_vec[0] = bf128_add(bf128_mul(q1, q2), bf128_mul(q3, bf_delta));

    for (unsigned int i = 1; i < n - 1; i++) {
        bf128_t q1 = bf_q[n + i - 1];
        bf128_t q2 = temp_q[i + 1];
        bf128_t q_or = bf_q[n + i];
        bf128_t q3 = bf128_add(q_or, bf128_add(q1, q2));

        b_vec[i] = bf128_add(bf128_mul(q1, q2), bf128_mul(q3, bf_delta));
    }

    free(temp_q);
}

void leave_prove_128(const uint8_t* w, const bf128_t* bf_v,
                     bf128_t* a_0_vec, bf128_t* a_1_vec, bf128_t* a_2_vec) {
    uint8_t* k = (uint8_t*)malloc((10 + 1) * 128 / 8);
    bf128_t* vk = (bf128_t*)faest_aligned_alloc(BF128_ALIGN, sizeof(bf128_t) * ((10 + 1) * 128));
    
    // w : t xor c || c || roundkeys || states
    aes_key_schedule_constraints_Mkey_0_128(
        w + 32, w + 16, bf_v + 256, bf_v + 128, a_0_vec, a_1_vec, a_2_vec, k, vk);
    
    uint8_t* in = (uint8_t*)malloc(16);
    uint8_t* out = (uint8_t*)malloc(16);
    bf128_t* bf_in = (bf128_t*)malloc(128 * sizeof(bf128_t));
    bf128_t* bf_out = (bf128_t*)malloc(128 * sizeof(bf128_t));
    
    xor_u8_array(w, w + 128 / 8, in, 16);
    xor_u8_array(w + AES_128_1 / 8, in, out, 16);
    
    for (int j = 0; j < 128; j++) {
        bf_in[j] = bf128_add(bf_v[j], bf_v[128 + j]);
        bf_out[j] = bf128_add(bf_v[AES_128_1 + j], bf_in[j]);
    }
    
    aes_enc_constraints_Mkey_0_128(in, out, w + (256 + 320) / 8,
                                   bf_v + 256 + 320, k, vk, bf_in, bf_out,
                                   a_0_vec + 80, a_1_vec + 80, a_2_vec + 80, 1, 1);
    
    faest_aligned_free(vk);
    free(k);
    free(in);
    free(out);
    free(bf_in);
    free(bf_out);
}

void leave_verify_128(const bf128_t* bf_q, const uint8_t* delta,
                      bf128_t* b_vec) {
    bf128_t* qk = (bf128_t*)faest_aligned_alloc(BF128_ALIGN, sizeof(bf128_t) * ((10 + 1) * 128));
    
    aes_key_schedule_constraints_Mkey_1_128(bf_q + 256, bf_q + 128, delta,
                                            b_vec, qk);
    
    bf128_t* bf_in = (bf128_t*)malloc(128 * sizeof(bf128_t));
    bf128_t* bf_out = (bf128_t*)malloc(128 * sizeof(bf128_t));
    
    for (int j = 0; j < 128; j++) {
        bf_in[j] = bf128_add(bf_q[j], bf_q[128 + j]);
        bf_out[j] = bf128_add(bf_q[AES_128_1 + j], bf_in[j]);
    }
    
    aes_enc_constraints_Mkey_1_128(NULL, NULL, bf_q + FAEST_128F_Lke + 128, qk,
                                   delta, bf_in, bf_out, b_vec + 80, 1, 1);
    
    faest_aligned_free(qk);
    free(bf_in);
    free(bf_out);
}

void merkle_tree_prove_128(const uint8_t* w, const bf128_t* bf_v,
                           bf128_t* a_0_vec, bf128_t* a_1_vec, bf128_t* a_2_vec,
                           unsigned int height) {
    uint8_t* k = (uint8_t*)malloc((10 + 1) * 128 / 8);
    bf128_t* vk = (bf128_t*)faest_aligned_alloc(BF128_ALIGN, sizeof(bf128_t) * ((10 + 1) * 128));
    uint8_t* origin_key = (uint8_t*)malloc(128 / 8);
    uint8_t* in = (uint8_t*)malloc(128 / 8);
    uint8_t* out = (uint8_t*)malloc(128 / 8);
    bf128_t* origin_key_v = (bf128_t*)malloc(128 * sizeof(bf128_t));
    bf128_t* bf_in = (bf128_t*)malloc(128 * sizeof(bf128_t));
    bf128_t* bf_out = (bf128_t*)malloc(128 * sizeof(bf128_t));
    
    for (int i = 0; i < height; i++) {
        xor_u8_array(w + MERKLE_TREE_128 * i / 8,
                     w + MERKLE_TREE_128 * i / 8 + 3 * 16UL, origin_key, 16UL);
        xor_u8_array(w + MERKLE_TREE_128 * i / 8 + 16UL,
                     w + 3 * 16UL + MERKLE_TREE_128 * i / 8, in, 16UL);
        xor_u8_array(w + MERKLE_TREE_128 / 8 + MERKLE_TREE_128 * i / 8, in, out,
                     16UL);
        
        for (int j = 0; j < 128; j++) {
            origin_key_v[j] =
                bf128_add(bf_v[j + MERKLE_TREE_128 * i],
                          bf_v[128 * 3 + j + MERKLE_TREE_128 * i]);
            bf_in[j] = bf128_add(bf_v[128 + j + MERKLE_TREE_128 * i],
                                 bf_v[128 * 3 + j + MERKLE_TREE_128 * i]);
            bf_out[j] = bf128_add(
                bf_v[MERKLE_TREE_128 + j + MERKLE_TREE_128 * i], bf_in[j]);
        }
        
        aes_key_schedule_constraints_Mkey_0_128(
            w + i * MERKLE_TREE_128 / 8 + 4 * 16, origin_key,
            bf_v + i * MERKLE_TREE_128 + 4 * 128, origin_key_v,
            a_0_vec + 400 * i, a_1_vec + 400 * i, a_2_vec + 400 * i, k, vk);
        
        aes_enc_constraints_Mkey_0_128(
            in, out, w + i * MERKLE_TREE_128 / 8 + 4 * 16 + 10 * 4,
            bf_v + i * MERKLE_TREE_128 + 4 * 128 + 10 * 32, k, vk, bf_in,
            bf_out, a_0_vec + 400 * i + 80, a_1_vec + 400 * i + 80, a_2_vec + 400 * i + 80, 1, 1);
    }
    
    faest_aligned_free(vk);
    free(k);
    free(origin_key);
    free(in);
    free(out);
    free(origin_key_v);
    free(bf_in);
    free(bf_out);
}

void merkle_tree_verify_128(const bf128_t* bf_q, const uint8_t* delta,
                            bf128_t* b_vec, unsigned int height) {
    bf128_t* qk = (bf128_t*)faest_aligned_alloc(BF128_ALIGN, sizeof(bf128_t) * ((10 + 1) * 128));
    bf128_t* origin_key_q = (bf128_t*)malloc(128 * sizeof(bf128_t));
    bf128_t* bf_in = (bf128_t*)malloc(128 * sizeof(bf128_t));
    bf128_t* bf_out = (bf128_t*)malloc(128 * sizeof(bf128_t));
    
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < 128; j++) {
            origin_key_q[j] =
                bf128_add(bf_q[j + MERKLE_TREE_128 * i],
                          bf_q[128 * 3 + j + MERKLE_TREE_128 * i]);
            bf_in[j] = bf128_add(bf_q[128 + j + MERKLE_TREE_128 * i],
                                 bf_q[128 * 3 + j + MERKLE_TREE_128 * i]);
            bf_out[j] = bf128_add(
                bf_q[MERKLE_TREE_128 + j + MERKLE_TREE_128 * i], bf_in[j]);
        }
        
        aes_key_schedule_constraints_Mkey_1_128(
            bf_q + i * MERKLE_TREE_128 + 4 * 128, origin_key_q, delta,
            b_vec + i * 400, qk);
        
        aes_enc_constraints_Mkey_1_128(
            NULL, NULL, bf_q + i * MERKLE_TREE_128 + 4 * 128 + 10 * 32, qk,
            delta, bf_in, bf_out, b_vec + 80 + i * 400, 1, 1);
    }
    
    free(origin_key_q);
    free(bf_in);
    free(bf_out);
    faest_aligned_free(qk);
}
