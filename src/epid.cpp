#include "epid.h"

#include <omp.h>

#include <chrono>
#include <cstdint>
#include <random>

#include "vole_parallel.h"
#include "witness_prove.hpp"

// Global fast random number generator - initialized once
static std::mt19937 global_rng(std::random_device{}());
static std::uniform_int_distribution<uint8_t> byte_dist(0, 255);

void hash_func(const std::vector<uint8_t>& m_1, const std::vector<uint8_t>& m_2,
               std::vector<uint8_t>& out) {
    EVP_CIPHER_CTX* ctx = EVP_CIPHER_CTX_new();
    EVP_EncryptInit_ex(ctx, EVP_aes_128_ecb(), NULL, m_1.data(), NULL);
    EVP_CIPHER_CTX_set_padding(ctx, 0);

    out.resize(16);
    int len;
    EVP_EncryptUpdate(ctx, out.data(), &len, m_2.data(), 16);
    EVP_EncryptFinal_ex(ctx, out.data() + len, &len);
    EVP_CIPHER_CTX_free(ctx);

    for (int i = 0; i < 16; i++) {
        out[i] ^= m_2[i];
    }
}

void hash_func_256(const std::vector<uint8_t>& m_1,
                   const std::vector<uint8_t>& m_2, std::vector<uint8_t>& out) {
    aes_round_keys_t round_keys;
    rijndael256_init_round_keys(&round_keys, m_1.data());

    out.resize(32);
    rijndael256_encrypt_block(&round_keys, m_2.data(), out.data());

    for (int i = 0; i < 32; i++) {
        out[i] ^= m_2[i];
    }
}

void gen_random(std::vector<uint8_t>& data) {
    for (auto& byte : data) {
        byte = byte_dist(global_rng);
    }
}

void epid::hash_mu(const std::vector<uint8_t>& msg, std::vector<uint8_t>& mu) {
    H1_context_t h1_ctx;
    H1_init(&h1_ctx, params_.lambda);
    H1_update(&h1_ctx, msg.data(), msg.size());
    H1_final(&h1_ctx, mu.data(), 2 * params_.lambda_bytes);
}

void epid::hash_challenge_1(const std::vector<uint8_t>& mu,
                            const std::vector<uint8_t>& hcom,
                            const std::vector<uint8_t>& c,
                            const std::vector<uint8_t>& iv,
                            std::vector<uint8_t>& chall_1, unsigned int ell,
                            unsigned int tau) {
    const unsigned int ell_hat_bytes =
        ell / 8 + params_.lambda_bytes * 3 + UNIVERSAL_HASH_B;
    H2_context_t h2_ctx;
    H2_init(&h2_ctx, params_.lambda);
    H2_update(&h2_ctx, mu.data(), params_.lambda_bytes * 2);
    H2_update(&h2_ctx, hcom.data(), params_.lambda_bytes * 2);
    H2_update(&h2_ctx, c.data(), ell_hat_bytes * (tau - 1));
    H2_update(&h2_ctx, iv.data(), IV_SIZE);
    H2_final(&h2_ctx, chall_1.data(), 5 * params_.lambda_bytes + 8);
}

void epid::hash_challenge_2(std::vector<uint8_t>& chall_2,
                            const std::vector<uint8_t>& chall_1,
                            const std::vector<uint8_t>& u_tilde,
                            const std::vector<uint8_t>& h_v,
                            const std::vector<uint8_t>& d, unsigned int lambda,
                            unsigned int ell) {
    const unsigned int lambda_bytes = lambda / 8;
    const unsigned int ell_bytes = ell / 8;
    const unsigned int u_tilde_bytes = lambda_bytes + UNIVERSAL_HASH_B;

    H2_context_t h2_ctx_1;
    H2_init(&h2_ctx_1, lambda);
    H2_update(&h2_ctx_1, chall_1.data(), 5 * lambda_bytes + 8);
    H2_update(&h2_ctx_1, u_tilde.data(), u_tilde_bytes);
    H2_update(&h2_ctx_1, h_v.data(), 2 * lambda_bytes);
    H2_update(&h2_ctx_1, d.data(), ell_bytes);
    H2_final(&h2_ctx_1, chall_2.data(), 3 * lambda_bytes + 8);
}

void epid::hash_challenge_3(std::vector<uint8_t>& chall_3,
                            const std::vector<uint8_t>& chall_2,
                            const std::vector<uint8_t>& a_0_tilde,
                            const std::vector<uint8_t>& a_1_tilde,
                            const std::vector<uint8_t>& a_2_tilde,
                            unsigned int lambda) {
    const unsigned int lambda_bytes = lambda / 8;

    H2_context_t h2_ctx_2;
    H2_init(&h2_ctx_2, lambda);
    H2_update(&h2_ctx_2, chall_2.data(), 3 * lambda_bytes + 8);
    H2_update(&h2_ctx_2, a_0_tilde.data(), lambda_bytes);
    H2_update(&h2_ctx_2, a_1_tilde.data(), lambda_bytes);
    H2_update(&h2_ctx_2, a_2_tilde.data(), lambda_bytes);
    H2_final(&h2_ctx_2, chall_3.data(), lambda_bytes);
}

bf128_t epid::zk_hash(const std::vector<uint8_t>& sd,
                      const std::vector<bf128_t>& x_0, const bf128_t& x_1) {
    bf128_t r_0 = bf128_load(sd.data());
    bf128_t r_1 = bf128_load(sd.data() + params_.lambda_bytes);
    bf128_t s = bf128_load(sd.data() + params_.lambda_bytes * 2);

    uint64_t tmp;
    memcpy(&tmp, sd.data() + 3 * params_.lambda_bytes, sizeof(uint64_t));
    bf128_t t = bf128_from_bf64(tmp);

    bf128_t h_0 = bf128_zero();
    bf128_t h_1 = bf128_zero();

    // Use simple serial computation - precomputation was too expensive
    bf128_t s_muti = s;
    bf128_t t_muti = t;

    for (size_t i = 0; i < x_0.size(); ++i) {
        h_0 = bf128_add(h_0, bf128_mul(s_muti, x_0[i]));
        h_1 = bf128_add(h_1, bf128_mul(t_muti, x_0[i]));
        s_muti = bf128_mul(s_muti, s);
        t_muti = bf128_mul(t_muti, t);
    }

    bf128_t result = bf128_mul(r_0, h_0);
    result = bf128_add(result, bf128_mul(r_1, h_1));
    result = bf128_add(result, x_1);

    return result;
}

void epid::gen_rootkey_iv(const std::vector<uint8_t>& mu, unsigned int index,
                          std::vector<uint8_t>& rootkey,
                          std::vector<uint8_t>& iv) {
    H3_context_t h3_ctx;
    H3_init(&h3_ctx, params_.lambda);
    H3_update(&h3_ctx, skey_set[index].data(), params_.lambda_bytes);
    H3_update(&h3_ctx, mu.data(), params_.lambda_bytes * 2);
    H3_final(&h3_ctx, rootkey.data(), params_.lambda_bytes, iv.data());
}

int hash_1(const uint8_t* key, const uint8_t* input, uint8_t* output) {
    aes_round_keys_t round_keys;
    int ret = aes128_init_round_keys(&round_keys, key);
    ret |= aes128_encrypt_block(&round_keys, input, output);
    xor_u8_array(input, output, output, 128 / 8);
    return ret == 0;
}

int hash_1_256(const uint8_t* key, const uint8_t* input, uint8_t* output) {
    aes_round_keys_t round_keys;
    int ret = rijndael256_init_round_keys(&round_keys, key);
    ret |= rijndael256_encrypt_block(&round_keys, input, output);
    xor_u8_array(input, output, output, 256 / 8);
    return ret == 0;
}

void epid::gen_sk() {
    std::vector<uint8_t> sk(lambda_bytes_);

    for (size_t i = 0; i < lambda_bytes_; ++i) {
        sk[i] = byte_dist(global_rng);
    }

    skey_set.push_back(std::move(sk));
}

void epid::gen_challenge() {
    std::vector<uint8_t> challenge(lambda_bytes_);

    for (size_t i = 0; i < lambda_bytes_; ++i) {
        challenge[i] = byte_dist(global_rng);
    }

    challenge_set.push_back(std::move(challenge));
}

void epid::gen_x() {
    std::vector<uint8_t> x;
    x.resize(lambda_bytes_);
    hash_func_256(challenge_set[challenge_set.size() - 1],
                  t_set[t_set.size() - 1], x);
    x_set.push_back(std::move(x));
}

void epid::cal_t(const std::vector<uint8_t>& sk,
                 const std::vector<uint8_t>& c_i) {
    std::vector<uint8_t> t(lambda_bytes_);
    hash_func_256(sk, c_i, t);
    t_set.push_back(std::move(t));
}

void epid::cal_t_256(const std::vector<uint8_t>& sk,
                     const std::vector<uint8_t>& c_i) {
    std::vector<uint8_t> t(lambda_bytes_);
    hash_func_256(sk, c_i, t);
    t_set.push_back(std::move(t));
}

void epid::gen_tree() {
    tree_ = std::vector<std::vector<uint8_t>>(
        2 * member_num_, std::vector<uint8_t>(lambda_bytes_));

#pragma omp parallel for schedule(static)
    for (int i = member_num_; i < 2 * member_num_; i++) {
        tree_[i] = x_set[i - member_num_];
    }

    for (int level = 1; level < log2(member_num_) + 1; level++) {
        int level_start = member_num_ >> level;
        int level_end = member_num_ >> (level - 1);

#pragma omp parallel for schedule(static)
        for (int index = level_start; index < level_end; index++) {
            hash_func_256(tree_[2 * index], tree_[2 * index + 1], tree_[index]);
        }
    }
}

void epid::issue_join_parallel(int index) {
    thread_local std::mt19937 local_rng(std::random_device{}() + index);
    thread_local std::uniform_int_distribution<uint8_t> local_byte_dist(0, 255);

    std::vector<uint8_t> sk(lambda_bytes_);
    for (size_t i = 0; i < lambda_bytes_; ++i) {
        sk[i] = local_byte_dist(local_rng);
    }
    skey_set[index] = std::move(sk);

    std::vector<uint8_t> challenge(lambda_bytes_);
    for (size_t i = 0; i < lambda_bytes_; ++i) {
        challenge[i] = local_byte_dist(local_rng);
    }
    challenge_set[index] = std::move(challenge);

    std::vector<uint8_t> t(lambda_bytes_);
    hash_func_256(skey_set[index], challenge_set[index], t);
    t_set[index] = std::move(t);

    std::vector<uint8_t> x(lambda_bytes_);
    hash_func_256(challenge_set[index], t_set[index], x);
    x_set[index] = std::move(x);
}

void epid::issue_join() {
    gen_sk();

    gen_challenge();

    cal_t(skey_set[skey_set.size() - 1],
          challenge_set[challenge_set.size() - 1]);
    gen_x();
}

// Generate SRL (Signature Revocation List) with random revoked signatures
void epid::generate_srl(unsigned int srl_size) {
    // Generate a random secret key (different from existing ones)
    std::vector<uint8_t> random_sk(lambda_bytes_);
    bool sk_unique;
    do {
        sk_unique = true;
        for (auto& byte : random_sk) {
            byte = byte_dist(global_rng);
        }

        // Check if it's different from existing secret keys
        for (const auto& existing_sk : skey_set) {
            if (std::equal(random_sk.begin(), random_sk.end(),
                           existing_sk.begin())) {
                sk_unique = false;
                break;
            }
        }
    } while (!sk_unique);

    for (unsigned int i = 0; i < srl_size; i++) {
        signature_t revoked_sig;
        revoked_sig.r.resize(lambda_bytes_);
        for (auto& byte : revoked_sig.r) {
            byte = byte_dist(global_rng);
        }

        // Calculate t = f(sk, r) using AES
        revoked_sig.t.resize(lambda_bytes_);
        hash_func(random_sk, revoked_sig.r, revoked_sig.t);

        srl.push_back(revoked_sig);
    }

    std::cout << "Generated SRL with " << srl_size << " revoked signatures"
              << std::endl;
}

// Generate SRL (Signature Revocation List) with random revoked signatures for
// 256-bit
void epid::generate_srl_256(unsigned int srl_size) {
    // Generate a random secret key (different from existing ones)
    std::vector<uint8_t> random_sk(lambda_bytes_);
    bool sk_unique;
    do {
        sk_unique = true;
        for (auto& byte : random_sk) {
            byte = byte_dist(global_rng);
        }

        // Check if it's different from existing secret keys
        for (const auto& existing_sk : skey_set) {
            if (std::equal(random_sk.begin(), random_sk.end(),
                           existing_sk.begin())) {
                sk_unique = false;
                break;
            }
        }
    } while (!sk_unique);

    for (unsigned int i = 0; i < srl_size; i++) {
        signature_t revoked_sig;
        revoked_sig.r.resize(lambda_bytes_);
        for (auto& byte : revoked_sig.r) {
            byte = byte_dist(global_rng);
        }

        // Calculate t = f(sk, r) using Rijndael-256
        revoked_sig.t.resize(lambda_bytes_);
        hash_func_256(random_sk, revoked_sig.r, revoked_sig.t);

        srl.push_back(revoked_sig);
    }

    std::cout << "Generated SRL with " << srl_size
              << " revoked signatures (256-bit)" << std::endl;
}

void epid::gen_witness_merkle_tree(uint8_t* witness,
                                   unsigned int signer_index) {
    std::vector<uint8_t> tmp;
    tmp.resize(lambda_bytes_);
    unsigned int index = signer_index + member_num_;

    while (index != 1) {
        if (index % 2) {
            hash_1_witness_multiplexer(tree_[index], tree_[index - 1], s_1_,
                                       witness);
            witness += MERKLE_TREE_256 / 8;

            index = index / 2;
        } else {
            hash_1_witness_multiplexer(tree_[index], tree_[index + 1], s_0_,
                                       witness);
            witness += MERKLE_TREE_256 / 8;
            index = index / 2;
        }
    }
    memcpy(witness, tree_[index].data(), lambda_bytes_);
}

void epid::hash_1_witness_multiplexer(const std::vector<uint8_t>& input_0,
                                      const std::vector<uint8_t>& input_1,
                                      const std::vector<uint8_t>& s_byte,
                                      uint8_t* witness) {
    std::vector<uint8_t> key_bytes(lambda_bytes_, 0);
    std::vector<uint8_t> msg_bytes(lambda_bytes_, 0);
    memcpy(witness, input_0.data(), lambda_bytes_);
    memcpy(witness + lambda_bytes_, input_1.data(), lambda_bytes_);
    memcpy(witness + lambda_bytes_ * 2, s_byte.data(), lambda_bytes_);
    witness += lambda_bytes_ * 3;
    field::GF2_256 i_0, i_1, s, s_, key, msg;

    i_0.from_bytes(input_0.data());
    i_1.from_bytes(input_1.data());
    s.from_bytes(s_byte.data());
    s_ = (i_0 + i_1) * s;

    s_.to_bytes(witness);
    witness += lambda_bytes_;

    key = s_ + i_0;
    msg = s_ + i_1;
    key.to_bytes(key_bytes.data());
    msg.to_bytes(msg_bytes.data());
    aes_extend_witness(key_bytes.data(), msg_bytes.data(), witness, 0, 1, 0,
                       lambda_);
}

static void or_extend_witness(uint8_t* witness, uint8_t* d,
                              unsigned int lambda) {
    unsigned int lambda_bytes = lambda / 8;
    unsigned int d_bits = lambda_bytes * 8;

    uint8_t bit0 = ptr_get_bit(d, 0);
    uint8_t bit1 = ptr_get_bit(d, 1);
    uint8_t result = bit0 | bit1;
    setbit(&witness[0], 0, result);

    for (unsigned int i = 1; i < d_bits - 1; i++) {
        uint8_t prev_witness_bit = ptr_get_bit(witness, i - 1);
        uint8_t d_bit = ptr_get_bit(d, i + 1);
        uint8_t or_result = prev_witness_bit | d_bit;
        setbit(&witness[i / 8], i % 8, or_result);
    }
}

// witness for relation : for any sig_j \in SRL, t_j \neq f(sk_i,r_j)
//   sk_i || key_expand || aes_mid_state_1 || aes_out_1 || or_mid_state_1 ||
//   .... ||  aes_mid_state_j || aes_out_j || or_mid_state_j || ...
static void srl_extend_witness(std::vector<signature_t>& revoke_list,
                               std::vector<uint8_t>& sk, uint8_t* witness,
                               unsigned int lambda) {
// Parallel processing of SRL list - each entry works on independent memory
// regions
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < revoke_list.size(); i++) {
        const auto& sig = revoke_list[i];

        // Calculate the witness offset for this SRL entry
        uint8_t* current_witness = witness + i * ((14 * 256 + 256) / 8);

        // Local buffer for this thread
        std::vector<uint8_t> d(lambda / 8);

        aes_extend_witness(sk.data(), sig.r.data(), current_witness, 0, 0, 1,
                           lambda);
        xor_u8_array(current_witness + 14 * 256 / 8 - lambda / 8, sig.r.data(),
                     d.data(),
                     lambda / 8);  // rijndael out xor with input : f output
        xor_u8_array(d.data(), sig.t.data(), d.data(), lambda / 8);
        or_extend_witness(current_witness + 14 * 256 / 8, d.data(), lambda);
    }
}

void epid::epid_sign(const std::vector<uint8_t>& msg, signature_t* sig,
                     unsigned int signer_index) {
    const unsigned int ell =
        256 + 14 * 64 + 13 * 256 + 13 * 256 + 256 + 256 + 14 * 64 + 13 * 256 +
        (256 * 3 + 256 + 14 * 64 + 256 * 13) * log2(member_num_) + 256 +
        (14 * 256 + 256) * srl.size();
    const unsigned int multi_num =
        1120 + 896 + 1120 + 1120 * log2(member_num_) + (896 + 255) * srl.size();
    const unsigned int ell_bytes = ell / 8;
    const unsigned int ell_hat =
        ell + params_.lambda * 3 + UNIVERSAL_HASH_B_BITS;
    const unsigned int ell_hat_bytes = ell_hat / 8;
    auto start_time = std::chrono::high_resolution_clock::now();
    std::cout << "\n=== EPID Sign Performance Analysis ===" << std::endl;

    std::vector<uint8_t> mu(2 * params_.lambda_bytes);
    hash_mu(msg, mu);

    std::vector<uint8_t> rootkey(params_.lambda_bytes);
    sig->iv.resize(iv_size_);
    gen_rootkey_iv(mu, signer_index, rootkey, sig->iv);

    std::vector<uint8_t> hcom(params_.lambda_bytes * 2);
    std::vector<vec_com_t> vecCom(params_.tau);
    std::vector<uint8_t> u(ell_hat_bytes);
    std::vector<uint8_t*> V(params_.lambda);
    V[0] = new uint8_t[params_.lambda * ell_hat_bytes];
    for (unsigned int i = 1; i < params_.lambda; ++i) {
        V[i] = V[0] + i * ell_hat_bytes;
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    sig->c.resize((params_.tau - 1) * ell_hat_bytes);

    auto vole_detailed_start = std::chrono::high_resolution_clock::now();
    vole_commit_parallel(rootkey.data(), sig->iv.data(), ell_hat, &(params_),
                         hcom.data(), vecCom.data(), sig->c.data(), u.data(),
                         V.data());
    auto vole_detailed_end = std::chrono::high_resolution_clock::now();

    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "[Sign] vole_commit: "
              << std::chrono::duration<double, std::milli>(t2 - t1).count()
              << " ms" << std::endl;
    std::vector<uint8_t> chall_1(5 * params_.lambda_bytes + 8);
    hash_challenge_1(mu, hcom, sig->c, sig->iv, chall_1, ell, params_.tau);

    sig->u_tilde.resize(params_.lambda_bytes + UNIVERSAL_HASH_B);
    vole_hash(sig->u_tilde.data(), chall_1.data(), u.data(), ell,
              params_.lambda);

    t1 = std::chrono::high_resolution_clock::now();
    std::vector<uint8_t> h_v(params_.lambda_bytes * 2);
    {
        std::vector<std::vector<uint8_t>> V_tilde_results(params_.lambda);
        for (auto& v : V_tilde_results) {
            v.resize(params_.lambda_bytes + UNIVERSAL_HASH_B);
        }

        auto vole_start = std::chrono::high_resolution_clock::now();

#pragma omp parallel for schedule(static)
        for (unsigned int i = 0; i < params_.lambda; ++i) {
            vole_hash(V_tilde_results[i].data(), chall_1.data(), V[i], ell,
                      params_.lambda);
        }
        auto vole_end = std::chrono::high_resolution_clock::now();
        std::cout << "    [h_v] parallel vole_hash (" << params_.lambda
                  << " calls): "
                  << std::chrono::duration<double, std::milli>(vole_end -
                                                               vole_start)
                         .count()
                  << " ms" << std::endl;

        auto h1_start = std::chrono::high_resolution_clock::now();
        // 串行更新H1上下文
        H1_context_t h1_ctx_1;
        H1_init(&h1_ctx_1, params_.lambda);
        for (unsigned int i = 0; i < params_.lambda; ++i) {
            H1_update(&h1_ctx_1, V_tilde_results[i].data(),
                      params_.lambda_bytes + UNIVERSAL_HASH_B);
        }
        H1_final(&h1_ctx_1, h_v.data(), params_.lambda_bytes * 2);
        auto h1_end = std::chrono::high_resolution_clock::now();
        std::cout << "    [h_v] H1_update/final: "
                  << std::chrono::duration<double, std::milli>(h1_end -
                                                               h1_start)
                         .count()
                  << " ms" << std::endl;
    }
    t2 = std::chrono::high_resolution_clock::now();
    std::cout << "[Sign] h_v computation (" << params_.lambda
              << " vole_hash calls): "
              << std::chrono::duration<double, std::milli>(t2 - t1).count()
              << " ms" << std::endl;
    std::vector<uint8_t> witness(ell_bytes, 0);
    unsigned int index = 0;

    t1 = std::chrono::high_resolution_clock::now();
    // t = f(sk_i, r)
    sig->r.resize(lambda_bytes_);
    gen_random(sig->r);
    aes_extend_witness(skey_set[signer_index].data(), sig->r.data(),
                       witness.data(), 1, 1, 1, lambda_);
    index += 256 + 14 * 64 + 13 * 256;
    sig->t.resize(lambda_bytes_);
    memcpy(sig->t.data(), witness.data() + index / 8, lambda_bytes_);

    // t_i_join = f(sk_i, c_i)
    aes_extend_witness(skey_set[signer_index].data(),
                       challenge_set[signer_index].data(),
                       witness.data() + index / 8, 0, 0, 1, lambda_);
    index += 14 * 256;  // states + rijndael_out(t_i_join xor c_i)

    // x_i = f(c_i,t_i_join) // c_i
    aes_extend_witness(challenge_set[signer_index].data(),
                       t_set[signer_index].data(), witness.data() + index / 8,
                       1, 1, 0, lambda_);
    index += 256 + 14 * 64 + 13 * 256;  // c_i_join + roundkeys + states
    t2 = std::chrono::high_resolution_clock::now();
    std::cout << "[Sign] 3x aes_extend_witness calls: "
              << std::chrono::duration<double, std::milli>(t2 - t1).count()
              << " ms" << std::endl;

    t1 = std::chrono::high_resolution_clock::now();
    // x_i accumulated in accumulator
    gen_witness_merkle_tree(witness.data() + index / 8, signer_index);
    index += (256 * 3 + 256 + 14 * 64 + 256 * 13) * log2(member_num_) + 256;
    t2 = std::chrono::high_resolution_clock::now();
    std::cout << "[Sign] gen_witness_merkle_tree (height=" << log2(member_num_)
              << "): "
              << std::chrono::duration<double, std::milli>(t2 - t1).count()
              << " ms" << std::endl;

    t1 = std::chrono::high_resolution_clock::now();
    srl_extend_witness(srl, skey_set[signer_index], witness.data() + index / 8,
                       lambda_);
    t2 = std::chrono::high_resolution_clock::now();
    std::cout << "[Sign] srl_extend_witness (" << srl.size()
              << " SRL entries): "
              << std::chrono::duration<double, std::milli>(t2 - t1).count()
              << " ms" << std::endl;

    sig->d.resize(ell_bytes);
    xor_u8_array(witness.data(), u.data(), sig->d.data(), ell_bytes);

    std::vector<uint8_t> chall_2(3 * params_.lambda_bytes + 8);
    hash_challenge_2(chall_2, chall_1, sig->u_tilde, h_v, sig->d,
                     params_.lambda, ell);

    std::vector<uint8_t> A_0_tilde_bytes(params_.lambda_bytes);
    sig->A_1_tilde_bytes.resize(params_.lambda_bytes);
    sig->A_2_tilde_bytes.resize(params_.lambda_bytes);
    std::vector<bf128_t> a_0_vec;
    std::vector<bf128_t> a_1_vec;
    std::vector<bf128_t> a_2_vec;
    a_0_vec.resize(multi_num);
    a_1_vec.resize(multi_num);
    a_2_vec.resize(multi_num);
    std::vector<bf128_t> a_0_vec_test;
    std::vector<bf128_t> a_1_vec_test;
    a_0_vec_test.resize(multi_num);
    a_1_vec_test.resize(multi_num);

    t1 = std::chrono::high_resolution_clock::now();
    bf128_t* bf_v = column_to_row_major_and_shrink_V_128(V.data(), ell_hat);
    auto t1_5 = std::chrono::high_resolution_clock::now();
    std::cout << "[Sign] column_to_row_major_and_shrink_V_128: "
              << std::chrono::duration<double, std::milli>(t1_5 - t1).count()
              << " ms" << std::endl;

    witness_prove_256(witness.data(), bf_v, a_0_vec.data(), a_1_vec.data(),
                      a_2_vec.data(), chall_2.data(), lambda_, sig->t.data(),
                      sig->r.data(), log2(member_num_), ell_hat, srl);
    t2 = std::chrono::high_resolution_clock::now();
    std::cout << "[Sign] witness_prove_256: "
              << std::chrono::duration<double, std::milli>(t2 - t1_5).count()
              << " ms" << std::endl;

    // ZK hash operations
    auto zk_start = std::chrono::high_resolution_clock::now();
    bf128_t u_0_mask = bf128_load(u.data() + ell_bytes);
    bf128_t u_1_mask = bf128_load(u.data() + ell_bytes + params_.lambda_bytes);

    bf128_t a_0_tilde, a_1_tilde, a_2_tilde;

#pragma omp parallel sections num_threads(3)
    {
#pragma omp section
        {
            a_0_tilde = zk_hash(chall_2, a_0_vec, bf128_sum_poly(bf_v + ell));
        }
#pragma omp section
        {
            a_1_tilde =
                zk_hash(chall_2, a_1_vec,
                        bf128_add(u_0_mask,
                                  bf128_sum_poly(bf_v + ell + params_.lambda)));
        }
#pragma omp section
        {
            a_2_tilde = zk_hash(chall_2, a_2_vec, u_1_mask);
        }
    }
    auto zk_end = std::chrono::high_resolution_clock::now();
    std::cout
        << "[Sign] ZK hash operations (3 calls): "
        << std::chrono::duration<double, std::milli>(zk_end - zk_start).count()
        << " ms" << std::endl;

    bf128_store(A_0_tilde_bytes.data(), a_0_tilde);
    bf128_store(sig->A_1_tilde_bytes.data(), a_1_tilde);
    bf128_store(sig->A_2_tilde_bytes.data(), a_2_tilde);
    sig->chall_3.resize(params_.lambda_bytes);
    sig->pdec.resize(params_.tau);
    sig->com.resize(params_.tau);

    auto hash3_start = std::chrono::high_resolution_clock::now();
    hash_challenge_3(sig->chall_3, chall_2, A_0_tilde_bytes,
                     sig->A_1_tilde_bytes, sig->A_2_tilde_bytes,
                     params_.lambda);
    auto hash3_end = std::chrono::high_resolution_clock::now();
    std::cout << "[Sign] hash_challenge_3: "
              << std::chrono::duration<double, std::milli>(hash3_end -
                                                           hash3_start)
                     .count()
              << " ms" << std::endl;

    auto vector_open_start = std::chrono::high_resolution_clock::now();
    for (unsigned int i = 0; i < params_.tau; i++) {
        // Step 20
        uint8_t s_[12];
        ChalDec(sig->chall_3.data(), i, params_.k0, params_.tau0, params_.k1,
                params_.tau1, s_);

        // Step 21
        const unsigned int depth = i < params_.tau0 ? params_.k0 : params_.k1;
        sig->pdec[i] = new uint8_t[depth * params_.lambda_bytes];
        sig->com[i] = new uint8_t[2 * params_.lambda_bytes];
        vector_open(vecCom[i].k, vecCom[i].com, s_, sig->pdec[i], sig->com[i],
                    depth, params_.lambda_bytes);
        vec_com_clear(&vecCom[i]);
    }
    auto vector_open_end = std::chrono::high_resolution_clock::now();
    std::cout << "[Sign] vector_open operations (" << params_.tau
              << " iterations): "
              << std::chrono::duration<double, std::milli>(vector_open_end -
                                                           vector_open_start)
                     .count()
              << " ms" << std::endl;

    auto end_time = std::chrono::high_resolution_clock::now();
    std::cout << "=== EPID Sign Total Time: "
              << std::chrono::duration<double, std::milli>(end_time -
                                                           start_time)
                     .count()
              << " ms ===\n"
              << std::endl;
    delete[] V[0];
    faest_aligned_free(bf_v);
}

bool epid::epid_verify(const std::vector<uint8_t>& msg,
                       const signature_t* sig) {
    const unsigned int ell =
        256 + 14 * 64 + 13 * 256 + 13 * 256 + 256 + 256 + 14 * 64 + 13 * 256 +
        (256 * 3 + 256 + 14 * 64 + 256 * 13) * log2(member_num_) + 256 +
        (14 * 256 + 256) * srl.size();
    const unsigned int multi_num =
        1120 + 896 + 1120 + 1120 * log2(member_num_) + (896 + 255) * srl.size();
    const unsigned int ell_bytes = ell / 8;
    const unsigned int ell_hat =
        ell + params_.lambda * 3 + UNIVERSAL_HASH_B_BITS;
    const unsigned int ell_hat_bytes = ell_hat / 8;

    auto start_time = std::chrono::high_resolution_clock::now();
    std::cout << "\n=== EPID Verify Performance Analysis ===" << std::endl;

    std::vector<uint8_t> mu(2 * params_.lambda_bytes);
    hash_mu(msg, mu);

    std::vector<uint8_t*> Q(params_.lambda);
    std::vector<uint8_t> hcom(params_.lambda_bytes * 2);
    Q[0] = new uint8_t[params_.lambda * ell_hat_bytes];
    for (unsigned int i = 1; i < params_.lambda; ++i) {
        Q[i] = Q[0] + i * ell_hat_bytes;
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    vole_reconstruct_parallel(sig->iv.data(), sig->chall_3.data(),
                              sig->pdec.data(), sig->com.data(), hcom.data(),
                              Q.data(), ell_hat, &params_);
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "[Verify] vole_reconstruct: "
              << std::chrono::duration<double, std::milli>(t2 - t1).count()
              << " ms" << std::endl;

    std::vector<uint8_t> chall_1(5 * params_.lambda_bytes + 8);
    hash_challenge_1(mu, hcom, sig->c, sig->iv, chall_1, ell, params_.tau);

    std::vector<uint8_t*> Q_(params_.lambda);
    Q_[0] = new uint8_t[params_.lambda * ell_hat_bytes];
    for (unsigned int i = 1; i < params_.lambda; ++i) {
        Q_[i] = Q_[0] + i * ell_hat_bytes;
    }

    std::vector<uint8_t*> Dtilde(params_.lambda);
    Dtilde[0] =
        new uint8_t[params_.lambda * (params_.lambda_bytes + UNIVERSAL_HASH_B)];
    for (unsigned int i = 1; i < params_.lambda; ++i) {
        Dtilde[i] = Dtilde[0] + i * (params_.lambda_bytes + UNIVERSAL_HASH_B);
    }
    memset(Dtilde[0], 0,
           params_.lambda * (params_.lambda_bytes + UNIVERSAL_HASH_B));

    t1 = std::chrono::high_resolution_clock::now();
    unsigned int Dtilde_idx = 0;
    unsigned int q_idx = 0;
    for (unsigned int i = 0; i < params_.tau; i++) {
        const unsigned int depth = i < params_.tau0 ? params_.k0 : params_.k1;

        // Step 11
        uint8_t delta[12];
        ChalDec(sig->chall_3.data(), i, params_.k0, params_.tau0, params_.k1,
                params_.tau1, delta);
        // Step 16
        for (unsigned int j = 0; j != depth; ++j, ++Dtilde_idx) {
            // for scan-build
            assert(Dtilde_idx < params_.lambda);
            masked_xor_u8_array(Dtilde[Dtilde_idx], sig->u_tilde.data(),
                                Dtilde[Dtilde_idx], delta[j],
                                params_.lambda_bytes + UNIVERSAL_HASH_B);
        }

        if (i == 0) {
            // Step 8
            memcpy(Q_[q_idx], Q[q_idx], ell_hat_bytes * depth);
            q_idx += depth;
        } else {
            // Step 14
            for (unsigned int d = 0; d < depth; ++d, ++q_idx) {
                masked_xor_u8_array(Q[q_idx],
                                    sig->c.data() + (i - 1) * ell_hat_bytes,
                                    Q_[q_idx], delta[d], ell_hat_bytes);
            }
        }
    }
    t2 = std::chrono::high_resolution_clock::now();
    std::cout << "[Verify] Q matrix processing (" << params_.tau
              << " iterations): "
              << std::chrono::duration<double, std::milli>(t2 - t1).count()
              << " ms" << std::endl;

    t1 = std::chrono::high_resolution_clock::now();
    std::vector<uint8_t> h_v(params_.lambda_bytes * 2);
    {
        std::vector<std::vector<uint8_t>> Q_tilde_results(params_.lambda);
        for (auto& q : Q_tilde_results) {
            q.resize(params_.lambda_bytes + UNIVERSAL_HASH_B);
        }

        auto vole_start = std::chrono::high_resolution_clock::now();

#pragma omp parallel for schedule(static)
        for (unsigned int i = 0; i < params_.lambda; i++) {
            vole_hash(Q_tilde_results[i].data(), chall_1.data(), Q_[i], ell,
                      params_.lambda);
            xor_u8_array(Q_tilde_results[i].data(), Dtilde[i],
                         Q_tilde_results[i].data(),
                         params_.lambda_bytes + UNIVERSAL_HASH_B);
        }
        auto vole_end = std::chrono::high_resolution_clock::now();
        std::cout << "    [h_v] parallel vole_hash+xor (" << params_.lambda
                  << " calls): "
                  << std::chrono::duration<double, std::milli>(vole_end -
                                                               vole_start)
                         .count()
                  << " ms" << std::endl;

        auto h1_start = std::chrono::high_resolution_clock::now();

        H1_context_t h1_ctx_1;
        H1_init(&h1_ctx_1, params_.lambda);
        for (unsigned int i = 0; i < params_.lambda; i++) {
            H1_update(&h1_ctx_1, Q_tilde_results[i].data(),
                      params_.lambda_bytes + UNIVERSAL_HASH_B);
        }
        H1_final(&h1_ctx_1, h_v.data(), params_.lambda_bytes * 2);
        auto h1_end = std::chrono::high_resolution_clock::now();
        std::cout << "    [h_v] H1_update/final: "
                  << std::chrono::duration<double, std::milli>(h1_end -
                                                               h1_start)
                         .count()
                  << " ms" << std::endl;
    }
    delete[] Dtilde[0];
    t2 = std::chrono::high_resolution_clock::now();
    std::cout << "[Verify] h_v computation (" << params_.lambda
              << " vole_hash calls): "
              << std::chrono::duration<double, std::milli>(t2 - t1).count()
              << " ms" << std::endl;

    std::vector<uint8_t> chall_2(3 * params_.lambda_bytes + 8);
    hash_challenge_2(chall_2, chall_1, sig->u_tilde, h_v, sig->d,
                     params_.lambda, ell);

    for (unsigned int i = 0, col = 0; i < params_.tau; i++) {
        unsigned int depth = i < params_.tau0 ? params_.k0 : params_.k1;
        uint8_t decoded_challenge[MAX_DEPTH];
        ChalDec(sig->chall_3.data(), i, params_.k0, params_.tau0, params_.k1,
                params_.tau1, decoded_challenge);
        for (unsigned int j = 0; j < depth; j++, ++col) {
            if (decoded_challenge[j] == 1) {
                xor_u8_array(sig->d.data(), Q_[col], Q_[col], ell_bytes);
            }
        }
    }

    std::vector<uint8_t> A_0_tilde_bytes(params_.lambda_bytes);

    std::vector<bf128_t> b_vec(multi_num);
    std::vector<bf128_t> b_vec_test(multi_num);

    t1 = std::chrono::high_resolution_clock::now();
    bf128_t* bf_q = column_to_row_major_and_shrink_V_128(Q_.data(), ell_hat);
    auto t1_5 = std::chrono::high_resolution_clock::now();
    std::cout << "[Verify] column_to_row_major_and_shrink_V_128: "
              << std::chrono::duration<double, std::milli>(t1_5 - t1).count()
              << " ms" << std::endl;

    witness_verify_256(bf_q, sig->chall_3.data(), b_vec.data(), lambda_,
                       sig->t.data(), sig->r.data(), log2(member_num_), ell_hat,
                       srl);
    t2 = std::chrono::high_resolution_clock::now();
    std::cout << "[Verify] witness_verify_256: "
              << std::chrono::duration<double, std::milli>(t2 - t1_5).count()
              << " ms" << std::endl;

    auto zk_verify_start = std::chrono::high_resolution_clock::now();
    bf128_t bf_delta = bf128_load(sig->chall_3.data());

    bf128_t q_0 = bf128_sum_poly(bf_q + ell);
    bf128_t q_1 = bf128_sum_poly(bf_q + ell + params_.lambda);
    bf128_t b_tilde =
        zk_hash(chall_2, b_vec, bf128_add(q_0, bf128_mul(q_1, bf_delta)));

    bf128_t a_1_tilde = bf128_load(sig->A_1_tilde_bytes.data());
    bf128_t a_2_tilde = bf128_load(sig->A_2_tilde_bytes.data());
    bf128_t bf_delta_sq = bf128_mul(bf_delta, bf_delta);
    bf128_t a_0_tilde =
        bf128_add(bf128_mul(a_2_tilde, bf_delta_sq),
                  bf128_add(b_tilde, bf128_mul(a_1_tilde, bf_delta)));

    bf128_store(A_0_tilde_bytes.data(), a_0_tilde);
    auto zk_verify_end = std::chrono::high_resolution_clock::now();
    std::cout << "[Verify] ZK verification operations: "
              << std::chrono::duration<double, std::milli>(zk_verify_end -
                                                           zk_verify_start)
                     .count()
              << " ms" << std::endl;

    std::vector<uint8_t> chall_3(params_.lambda_bytes);

    auto final_hash_start = std::chrono::high_resolution_clock::now();
    hash_challenge_3(chall_3, chall_2, A_0_tilde_bytes, sig->A_1_tilde_bytes,
                     sig->A_2_tilde_bytes, params_.lambda);
    auto final_hash_end = std::chrono::high_resolution_clock::now();
    std::cout << "[Verify] final hash_challenge_3: "
              << std::chrono::duration<double, std::milli>(final_hash_end -
                                                           final_hash_start)
                     .count()
              << " ms" << std::endl;

    auto end_time = std::chrono::high_resolution_clock::now();
    std::cout << "=== EPID Verify Total Time: "
              << std::chrono::duration<double, std::milli>(end_time -
                                                           start_time)
                     .count()
              << " ms ===\n"
              << std::endl;

    delete[] Q[0];
    delete[] Q_[0];
    faest_aligned_free(bf_q);
    return memcmp(chall_3.data(), sig->chall_3.data(), params_.lambda_bytes) ==
           0;
}
