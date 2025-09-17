/*
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "aes_prove.h"
#include "parameters.h"
#include "universal_hashing.h"
#include "utils.h"
#include "vole.h"

static const bf8_t Rcon[30] = {
    0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36,
    0x6c, 0xd8, 0xab, 0x4d, 0x9a, 0x2f, 0x5e, 0xbc, 0x63, 0xc6,
    0x97, 0x35, 0x6a, 0xd4, 0xb3, 0x7d, 0xfa, 0xef, 0xc5, 0x91,
};

bf128_t* column_to_row_major_and_shrink_V_128(uint8_t** v, unsigned int ell) {
    // V is \hat \ell times \lambda matrix over F_2
    // v has \hat \ell rows, \lambda columns, storing in column-major order,
    // new_v has \ell + \lambda rows and \lambda columns storing in row-major
    // order
    bf128_t* new_v = faest_aligned_alloc(BF128_ALIGN, ell * sizeof(bf128_t));
    
    // Parallel processing of rows
    #pragma omp parallel for schedule(static)
    for (unsigned int row = 0; row < ell; ++row) {
        uint8_t new_row[BF128_NUM_BYTES] = {0};
        for (unsigned int column = 0; column != 128; ++column) {
            ptr_set_bit(new_row, ptr_get_bit(v[column], row), column);
        }
        new_v[row] = bf128_load(new_row);
    }

    return new_v;
}

// done
static void aes_key_schedule_forward_1(const uint8_t* x,
                                       const uint8_t* origin_key,
                                       uint8_t* out) {
    // Step: 1 skipped (sanity check)

    const unsigned int lambda = 128;
    const unsigned int R = 10;
    const unsigned int Nwd = 4;
    const unsigned int lambdaBytes = 16;

    const unsigned int out_len = (R + 1) * 128 / 8;
    // Step 3
    memcpy(out, origin_key, lambdaBytes);
    memset(out + lambdaBytes, 0, out_len - lambdaBytes);
    // Step: 4
    unsigned int i_wd = 0;
    // Step: 5..10
    for (unsigned int j = Nwd; j < 4 * (R + 1); j++) {
        if ((j % Nwd) == 0 || (Nwd > 6 && (j % Nwd) == 4)) {
            memcpy(out + 32 * j / 8, x + i_wd / 8, 4);
            i_wd += 32;
        } else {
            for (unsigned int i = 0; i < 32; i += 8) {
                // bit spliced
                out[(32 * j + i) / 8] |=
                    out[(32 * (j - Nwd) + i) / 8] ^ out[(32 * (j - 1) + i) / 8];
            }
        }
    }
}
// done
static void rijnd_key_schedule_forward_1(const uint8_t* x,
                                         const uint8_t* origin_key,
                                         uint8_t* out) {
    // Step: 1 skipped (sanity check)

    const unsigned int lambda = 256;
    const unsigned int R = 14;
    const unsigned int Nwd = 8;
    const unsigned int lambdaBytes = lambda / 8;

    const unsigned int out_len = (R + 1) * 256 / 8;
    // Step 3
    memcpy(out, origin_key, lambdaBytes);
    memset(out + lambdaBytes, 0, out_len - lambdaBytes);

    // Step: 4
    unsigned int i_wd = 0;
    // Step: 5..10
    for (unsigned int j = Nwd; j < Nwd * (R + 1); j++) {
        if ((j % Nwd) == 0 || (Nwd > 6 && (j % Nwd) == 4)) {
            memcpy(out + 32 * j / 8, x + i_wd / 8, 4);
            i_wd += 32;
        } else {
            for (unsigned int i = 0; i < 32; i += 8) {
                // bit spliced
                out[(32 * j + i) / 8] |=
                    out[(32 * (j - Nwd) + i) / 8] ^ out[(32 * (j - 1) + i) / 8];
            }
        }
    }
}

// done
static void aes_key_schedule_backward_1(const uint8_t* x, const uint8_t* xk,
                                        uint8_t* out) {
    // Step: 1 skipped (sanity check)

    const unsigned int lambda = 128;
    const unsigned int Ske = 40;

    // Step: 2
    unsigned int iwd = 0;
    unsigned int c = 0;
    bool rmvRcon = true;
    unsigned int ircon = 0;

    for (unsigned int j = 0; j < Ske; j++) {
        uint8_t x_tilde = x[j] ^ xk[iwd + c];
        if (/* Mtag == 0 && */ rmvRcon == true && c == 0) {
            // Steps 12 and 13, bitsliced; delta is always 0
            x_tilde ^= Rcon[ircon];
            ++ircon;
        }
        const uint8_t y_tilde =
            rotr8(x_tilde, 7) ^ rotr8(x_tilde, 5) ^ rotr8(x_tilde, 2);
        out[j] = y_tilde ^ 0x5;
        ++c;
        if (c == 4) {
            c = 0;
            iwd += 128 / 8;
        }
    }
}

// done
static void rijnd_key_schedule_backward_1(const uint8_t* x, const uint8_t* xk,
                                          uint8_t* out) {
    // Step: 1 skipped (sanity check)

    const unsigned int lambda = 256;
    const unsigned int Ske = 112;

    // Step: 2
    unsigned int iwd = 0;
    unsigned int c = 0;
    bool rmvRcon = true;
    unsigned int ircon = 0;

    for (unsigned int j = 0; j < Ske; j++) {
        uint8_t x_tilde = x[j] ^ xk[iwd + c];

        if (rmvRcon == true && c == 0) {
            x_tilde ^= Rcon[ircon];
            ++ircon;
        }

        const uint8_t y_tilde =
            rotr8(x_tilde, 7) ^ rotr8(x_tilde, 5) ^ rotr8(x_tilde, 2);
        out[j] = y_tilde ^ 0x5;

        ++c;
        if (c == 4) {
            c = 0;
            if (lambda == 192) {
                iwd += 192 / 8;
            } else {
                iwd += 128 / 8;
                if (lambda == 256) {
                    rmvRcon = !rmvRcon;
                }
            }
        }
    }
}

// lambda == 128 implementation
// done
static void aes_key_schedule_forward_128(const bf128_t* v,
                                         const bf128_t* origin_key_v,
                                         bf128_t* bf_out) {
    // Step: 1 sanity check (skipped)

    memcpy(bf_out, origin_key_v, 128 * sizeof(bf128_t));

    // Step: 4
    unsigned int i_wd = 0;
    // Step: 5..10
    for (unsigned int j = FAEST_128F_Nwd; j < 4 * (10 + 1); j++) {
        if ((j % FAEST_128F_Nwd) == 0 ||
            (FAEST_128F_Nwd > 6 && (j % FAEST_128F_Nwd) == 4)) {
            // copy all at once
            memcpy(bf_out + j * 32, v + i_wd, sizeof(bf128_t) * 32);
            i_wd += 32;
        } else {
            for (unsigned int i = 0; i < 32; i++) {
                bf_out[(32 * j) + i] =
                    bf128_add(bf_out[32 * (j - FAEST_128F_Nwd) + i],
                              bf_out[32 * (j - 1) + i]);
            }
        }
    }
}
// done
static void aes_key_schedule_backward_128(const bf128_t* v, const bf128_t* Vk,
                                          uint8_t Mtag, uint8_t Mkey,
                                          const uint8_t* delta,
                                          bf128_t* bf_out) {
    // Step: 1
    assert(!((Mtag == 1 && Mkey == 1) || (Mkey == 1 && delta == NULL)));

    const bf128_t bf_delta = delta ? bf128_load(delta) : bf128_zero();

    unsigned int iwd = 0;
    unsigned int c = 0;
    unsigned int ircon = 0;

    bf128_t bf_minus_mkey = bf128_from_bit(1 ^ Mkey);
    uint8_t minus_mtag = 1 ^ Mtag;
    bf128_t bf_mkey_times_delta = bf128_mul_bit(bf_delta, Mkey);
    bf_mkey_times_delta = bf128_add(bf_mkey_times_delta, bf_minus_mkey);

    for (unsigned int j = 0; j < FAEST_128F_Ske; j++) {
        // Step 7
        bf128_t bf_x_tilde[8];
        for (unsigned int i = 0; i < 8; i++) {
            bf_x_tilde[i] = bf128_add(v[8 * j + i], Vk[iwd + 8 * c + i]);
        }

        if (Mtag == 0 && c == 0) {
            // Step 9
            uint8_t r = Rcon[ircon];
            ircon = ircon + 1;

            bf128_t bf_r[8];
            for (unsigned int i = 0; i < 8; i++) {
                // Step 12
                bf_r[i] = bf128_mul_bit(bf_mkey_times_delta, get_bit(r, i));
                // Step 13
                bf_x_tilde[i] = bf128_add(bf_x_tilde[i], bf_r[i]);
            }
        }

        for (unsigned int i = 0; i < 8; ++i) {
            bf_out[i + 8 * j] = bf128_add(
                bf128_add(bf_x_tilde[(i + 7) % 8], bf_x_tilde[(i + 5) % 8]),
                bf_x_tilde[(i + 2) % 8]);
        }
        bf_out[0 + 8 * j] = bf128_add(
            bf_out[0 + 8 * j], bf128_mul_bit(bf_mkey_times_delta, minus_mtag));
        bf_out[2 + 8 * j] = bf128_add(
            bf_out[2 + 8 * j], bf128_mul_bit(bf_mkey_times_delta, minus_mtag));
        c = c + 1;

        if (c == 4) {
            c = 0;
            iwd += 128;
        }
    }
}
// done
void aes_key_schedule_constraints_Mkey_0_128(
    const uint8_t* w, const uint8_t* origin_key, const bf128_t* v,
    const bf128_t* origin_key_v, bf128_t* a_0_vec, bf128_t* a_1_vec,
    bf128_t* a_2_vec, uint8_t* k, bf128_t* vk) {
    // Step: 2
    aes_key_schedule_forward_1(w, origin_key, k);

    // Step: 3
    aes_key_schedule_forward_128(v, origin_key_v, vk);

    // Step: 4
    uint8_t w_dash[FAEST_128F_Ske];
    aes_key_schedule_backward_1(w, k, w_dash);

    // Step: 5
    bf128_t v_w_dash[FAEST_128F_Ske * 8];
    aes_key_schedule_backward_128(v, vk, 1, 0, NULL, v_w_dash);

    // Step: 6..8
    unsigned int iwd = 32 * (FAEST_128F_Nwd - 1);
    for (unsigned int j = 0; j < FAEST_128F_Ske / 4; j++) {
        bf128_t bf_k_hat[4];
        bf128_t bf_v_k_hat[4];
        bf128_t bf_w_dash_hat[4];
        bf128_t bf_v_w_dash_hat[4];
        for (unsigned int r = 0; r <= 3; r++) {
            // Step: 10..11
            bf_k_hat[(r + 3) % 4] =
                bf128_byte_combine_bits(k[(iwd + 8 * r) / 8]);
            bf_v_k_hat[(r + 3) % 4] = bf128_byte_combine(vk + (iwd + 8 * r));
            bf_w_dash_hat[r] =
                bf128_byte_combine_bits(w_dash[(32 * j + 8 * r) / 8]);
            bf_v_w_dash_hat[r] =
                bf128_byte_combine(v_w_dash + (32 * j + 8 * r));
        }
        // Step: 13..17
        for (unsigned int r = 0; r <= 3; r++) {
            unsigned int index = 4 * j + r;
            bf128_t v1 = bf_v_k_hat[r];
            bf128_t v2 = bf_v_w_dash_hat[r];
            bf128_t u1 = bf_k_hat[r];
            bf128_t u2 = bf_w_dash_hat[r];

            // a_0_vec[4 * j + r] = v1^2*v2
            bf128_t v1_squared = bf128_mul(v1, v1);
            a_0_vec[2 * index] = bf128_mul(v1_squared, v2);
            // a_1_vec[4 * j + r] = v1^2*u2
            a_1_vec[2 * index] = bf128_mul(v1_squared, u2);
            // a_2_vec[4 * j + r] = u1^2*v2 + v1
            bf128_t u1_squared = bf128_mul(u1, u1);
            a_2_vec[2 * index] = bf128_add(bf128_mul(u1_squared, v2), v1);

            // a_0_vec[4 * j + r] = v2^2*v1
            bf128_t v2_squared = bf128_mul(v2, v2);
            a_0_vec[2 * index + 1] = bf128_mul(v2_squared, v1);
            // a_1_vec[4 * j + r] = v2^2*u1
            a_1_vec[2 * index + 1] = bf128_mul(v2_squared, u1);
            // a_2_vec[4 * j + r] = u2^2*v1 + v2
            bf128_t u2_squared = bf128_mul(u2, u2);
            a_2_vec[2 * index + 1] = bf128_add(bf128_mul(u2_squared, v1), v2);
        }
        iwd = iwd + 128;
    }
}
// done
void aes_key_schedule_constraints_Mkey_1_128(const bf128_t* q,
                                             const bf128_t* origin_key_q,
                                             const uint8_t* delta,
                                             bf128_t* b_vec, bf128_t* qk) {
    // Step: 19..20
    aes_key_schedule_forward_128(q, origin_key_q, qk);
    bf128_t q_w_dash[FAEST_128F_Ske * 8];
    aes_key_schedule_backward_128(q, qk, 0, 1, delta, q_w_dash);

    const bf128_t bf_delta = bf128_load(delta);
    const bf128_t delta_squared = bf128_mul(bf_delta, bf_delta);

    // Step 23..24
    unsigned int iwd = 32 * (FAEST_128F_Nwd - 1);
    for (unsigned int j = 0; j < FAEST_128F_Ske / 4; j++) {
        bf128_t bf_q_hat_k[4];
        bf128_t bf_q_hat_w_dash[4];
        for (unsigned int r = 0; r <= 3; r++) {
            // Step: 25..26
            bf_q_hat_k[(r + 3) % 4] = bf128_byte_combine(qk + ((iwd + 8 * r)));
            bf_q_hat_w_dash[r] =
                bf128_byte_combine(q_w_dash + ((32 * j + 8 * r)));
        }
        for (unsigned int r = 0; r <= 3; r++) {
            unsigned int index = j * 4 + r;
            bf128_t q1 = bf_q_hat_k[r];
            bf128_t q2 = bf_q_hat_w_dash[r];

            // b_vec[2 * index] = q1^2*q2 - q1*delta^2
            bf128_t q1_squared = bf128_mul(q1, q1);
            b_vec[2 * index] = bf128_add(bf128_mul(q1_squared, q2),
                                         bf128_mul(q1, delta_squared));
            // b_vec[2 * index + 1] = q2^2*q1 - q2*delta^2
            bf128_t q2_squared = bf128_mul(q2, q2);
            b_vec[2 * index + 1] = bf128_add(bf128_mul(q2_squared, q1),
                                             bf128_mul(q2, delta_squared));
        }
        iwd = iwd + 128;
    }
}

static void rijnd_key_schedule_forward_256(const bf128_t* v,
                                           const bf128_t* origin_key_v,
                                           bf128_t* bf_out) {
    // Step: 1 sanity check (skipped)

    memcpy(bf_out, origin_key_v, sizeof(bf128_t) * FAEST_256F_LAMBDA);

    // Step: 4
    unsigned int i_wd = 0;
    // Step: 5..10
    for (unsigned int j = FAEST_256F_Nwd; j < 8 * (FAEST_256F_R + 1); j++) {
        if ((j % FAEST_256F_Nwd) == 0 ||
            (FAEST_256F_Nwd > 6 && (j % FAEST_256F_Nwd) == 4)) {
            memcpy(bf_out + j * 32, v + i_wd, sizeof(bf128_t) * 32);
            i_wd += 32;
        } else {
            for (unsigned int i = 0; i < 32; i++) {
                bf_out[(32 * j) + i] =
                    bf128_add(bf_out[32 * (j - FAEST_256F_Nwd) + i],
                              bf_out[32 * (j - 1) + i]);
            }
        }
    }
}

static void rijnd_key_schedule_backward_256(const bf128_t* v, const bf128_t* Vk,
                                            uint8_t Mtag, uint8_t Mkey,
                                            const uint8_t* delta,
                                            bf128_t* bf_out) {
    // Step: 1
    assert(!((Mtag == 1 && Mkey == 1) || (Mkey == 1 && delta == NULL)));

    unsigned int Ske = 112;
    unsigned int iwd = 0;
    unsigned int c = 0;
    bool rmvRcon = true;
    unsigned int ircon = 0;

    const bf128_t bf_delta = delta ? bf128_load(delta) : bf128_zero();
    const bf128_t bf_minus_mkey = bf128_from_bit(1 ^ Mkey);
    const uint8_t minus_mtag = 1 ^ Mtag;
    bf128_t bf_mkey_times_delta = bf128_mul_bit(bf_delta, Mkey);
    bf_mkey_times_delta = bf128_add(bf_mkey_times_delta, bf_minus_mkey);

    for (unsigned int j = 0; j < Ske; j++) {
        // Step 7
        bf128_t bf_x_tilde[8];
        for (unsigned int i = 0; i < 8; i++) {
            bf_x_tilde[i] = bf128_add(v[8 * j + i], Vk[iwd + 8 * c + i]);
        }

        if (Mtag == 0 && rmvRcon == true && c == 0) {
            // Step 9
            uint8_t r = Rcon[ircon];
            ircon = ircon + 1;

            bf128_t bf_r[8];
            for (unsigned int i = 0; i < 8; i++) {
                // Step 12
                bf_r[i] = bf128_mul_bit(bf_mkey_times_delta, get_bit(r, i));
                // Step 13
                bf_x_tilde[i] = bf128_add(bf_x_tilde[i], bf_r[i]);
            }
        }

        for (unsigned int i = 0; i < 8; ++i) {
            bf_out[i + 8 * j] = bf128_add(
                bf128_add(bf_x_tilde[(i + 7) % 8], bf_x_tilde[(i + 5) % 8]),
                bf_x_tilde[(i + 2) % 8]);
        }
        bf_out[0 + 8 * j] = bf128_add(
            bf_out[0 + 8 * j], bf128_mul_bit(bf_mkey_times_delta, minus_mtag));
        bf_out[2 + 8 * j] = bf128_add(
            bf_out[2 + 8 * j], bf128_mul_bit(bf_mkey_times_delta, minus_mtag));
        c = c + 1;

        if (c == 4) {
            c = 0;
            iwd += 128;
            rmvRcon = !rmvRcon;
        }
    }
}

void rijnd_key_schedule_constraints_Mkey_0_256(
    const uint8_t* w, const uint8_t* origin_key, const bf128_t* v,
    const bf128_t* origin_key_v, bf128_t* a_0_vec, bf128_t* a_1_vec,
    bf128_t* a_2_vec, uint8_t* k, bf128_t* vk) {
    // Step: 2
    rijnd_key_schedule_forward_1(w, origin_key, k);

    // Step: 3
    rijnd_key_schedule_forward_256(v, origin_key_v, vk);

    // Step: 4
    uint8_t w_dash[112];
    rijnd_key_schedule_backward_1(w, k, w_dash);

    // Step: 5
    bf128_t v_w_dash[112 * 8];
    rijnd_key_schedule_backward_256(v, vk, 1, 0, NULL, v_w_dash);

    // Step: 6..8
    bool rotate_word = true;
    unsigned int iwd = 32 * (FAEST_256F_Nwd - 1);
    for (unsigned int j = 0; j < 112 / 4; j++) {
        bf128_t bf_k_hat[4];
        bf128_t bf_v_k_hat[4];
        bf128_t bf_w_dash_hat[4];
        bf128_t bf_v_w_dash_hat[4];
        for (unsigned int r = 0; r <= 3; r++) {
            // Step: 10..11
            if (rotate_word) {
                bf_k_hat[(r + 3) % 4] =
                    bf128_byte_combine_bits(k[(iwd + 8 * r) / 8]);
                bf_v_k_hat[(r + 3) % 4] =
                    bf128_byte_combine(vk + (iwd + 8 * r));
                bf_w_dash_hat[r] =
                    bf128_byte_combine_bits(w_dash[(32 * j + 8 * r) / 8]);
                bf_v_w_dash_hat[r] =
                    bf128_byte_combine(v_w_dash + (32 * j + 8 * r));
            } else {
                bf_k_hat[r] = bf128_byte_combine_bits(k[(iwd + 8 * r) / 8]);
                bf_v_k_hat[r] = bf128_byte_combine(vk + (iwd + 8 * r));
                bf_w_dash_hat[r] =
                    bf128_byte_combine_bits(w_dash[(32 * j + 8 * r) / 8]);
                bf_v_w_dash_hat[r] =
                    bf128_byte_combine(v_w_dash + (32 * j + 8 * r));
            }
        }
        // Step: 13..17
        for (unsigned int r = 0; r <= 3; r++) {
            unsigned int index = 4 * j + r;
            bf128_t v1 = bf_v_k_hat[r];
            bf128_t v2 = bf_v_w_dash_hat[r];
            bf128_t u1 = bf_k_hat[r];
            bf128_t u2 = bf_w_dash_hat[r];

            // a_0_vec[4 * j + r] = v1^2*v2
            bf128_t v1_squared = bf128_mul(v1, v1);
            a_0_vec[2 * index] = bf128_mul(v1_squared, v2);
            // a_1_vec[4 * j + r] = v1^2*u2
            a_1_vec[2 * index] = bf128_mul(v1_squared, u2);
            // a_2_vec[4 * j + r] = u1^2*v2 + v1
            bf128_t u1_squared = bf128_mul(u1, u1);
            a_2_vec[2 * index] = bf128_add(bf128_mul(u1_squared, v2), v1);

            // a_0_vec[4 * j + r] = v2^2*v1
            bf128_t v2_squared = bf128_mul(v2, v2);
            a_0_vec[2 * index + 1] = bf128_mul(v2_squared, v1);
            // a_1_vec[4 * j + r] = v2^2*u1
            a_1_vec[2 * index + 1] = bf128_mul(v2_squared, u1);
            // a_2_vec[4 * j + r] = u2^2*v1 + v2
            bf128_t u2_squared = bf128_mul(u2, u2);
            a_2_vec[2 * index + 1] = bf128_add(bf128_mul(u2_squared, v1), v2);
        }
        iwd = iwd + 128;
        rotate_word = !rotate_word;
    }
}

void rijnd_key_schedule_constraints_Mkey_1_256(const bf128_t* q,
                                               const bf128_t* origin_key_q,
                                               const uint8_t* delta,
                                               bf128_t* b_vec, bf128_t* qk) {
    // Step: 19..20
    rijnd_key_schedule_forward_256(q, origin_key_q, qk);
    bf128_t q_w_dash[112 * 8];
    rijnd_key_schedule_backward_256(q, qk, 0, 1, delta, q_w_dash);

    const bf128_t bf_delta = bf128_load(delta);
    const bf128_t bf_delta_squared = bf128_mul(bf_delta, bf_delta);
    bool rotate_word = true;

    // Step 23..24
    unsigned int iwd = 32 * (FAEST_256F_Nwd - 1);
    for (unsigned int j = 0; j < 112 / 4; j++) {
        bf128_t bf_q_hat_k[4];
        bf128_t bf_q_hat_w_dash[4];
        for (unsigned int r = 0; r <= 3; r++) {
            // Step: 25..26
            if (rotate_word) {
                bf_q_hat_k[(r + 3) % 4] =
                    bf128_byte_combine(qk + ((iwd + 8 * r)));
                bf_q_hat_w_dash[r] =
                    bf128_byte_combine(q_w_dash + ((32 * j + 8 * r)));
            } else {
                bf_q_hat_k[r] = bf128_byte_combine(qk + ((iwd + 8 * r)));
                bf_q_hat_w_dash[r] =
                    bf128_byte_combine(q_w_dash + ((32 * j + 8 * r)));
            }
        }
        // Step: 27
        for (unsigned int r = 0; r <= 3; r++) {
            unsigned int index = j * 4 + r;
            bf128_t q1 = bf_q_hat_k[r];
            bf128_t q2 = bf_q_hat_w_dash[r];

            // b_vec[2 * index] = q1^2*q2 - q1*delta^2
            bf128_t q1_squared = bf128_mul(q1, q1);
            b_vec[2 * index] = bf128_add(bf128_mul(q1_squared, q2),
                                         bf128_mul(q1, bf_delta_squared));
            // b_vec[2 * index + 1] = q2^2*q1 - q2*delta^2
            bf128_t q2_squared = bf128_mul(q2, q2);
            b_vec[2 * index + 1] = bf128_add(bf128_mul(q2_squared, q1),
                                             bf128_mul(q2, bf_delta_squared));
        }
        iwd = iwd + 128;
        rotate_word = !rotate_word;
    }
}

// done
static void aes_enc_forward_128_1(const uint8_t* x, const uint8_t* xk,
                                  const uint8_t* in, bf128_t* bf_y) {
    for (unsigned int i = 0; i < 16; i++) {
        const uint8_t xin = in[i];
        bf_y[i] = bf128_add(bf128_byte_combine_bits(xin),
                            bf128_byte_combine_bits(xk[i]));
    }

    const bf128_t bf_two = bf128_byte_combine_bits(2);
    const bf128_t bf_three = bf128_byte_combine_bits(3);

    for (unsigned int j = 1; j < 10; j++) {
        for (unsigned int c = 0; c <= 3; c++) {
            const unsigned int ix = 128 * (j - 1) + 32 * c;
            const unsigned int ik = 128 * j + 32 * c;
            const unsigned int iy = 16 * j + 4 * c;

            bf128_t bf_x_hat[4];
            bf128_t bf_xk_hat[4];
            for (unsigned int r = 0; r <= 3; r++) {
                // Step: 12..13
                bf_x_hat[r] = bf128_byte_combine_bits(x[(ix + 8 * r) / 8]);
                bf_xk_hat[r] = bf128_byte_combine_bits(xk[(ik + 8 * r) / 8]);
            }

            // Step : 14
            bf_y[iy + 0] =
                bf128_add(bf_xk_hat[0], bf128_mul(bf_x_hat[0], bf_two));
            bf_y[iy + 0] =
                bf128_add(bf_y[iy + 0], bf128_mul(bf_x_hat[1], bf_three));
            bf_y[iy + 0] = bf128_add(bf_y[iy + 0], bf_x_hat[2]);
            bf_y[iy + 0] = bf128_add(bf_y[iy + 0], bf_x_hat[3]);

            // Step: 15
            bf_y[iy + 1] = bf128_add(bf_xk_hat[1], bf_x_hat[0]);
            bf_y[iy + 1] =
                bf128_add(bf_y[iy + 1], bf128_mul(bf_x_hat[1], bf_two));
            bf_y[iy + 1] =
                bf128_add(bf_y[iy + 1], bf128_mul(bf_x_hat[2], bf_three));
            bf_y[iy + 1] = bf128_add(bf_y[iy + 1], bf_x_hat[3]);

            // Step: 16
            bf_y[iy + 2] = bf128_add(bf_xk_hat[2], bf_x_hat[0]);
            bf_y[iy + 2] = bf128_add(bf_y[iy + 2], bf_x_hat[1]);
            bf_y[iy + 2] =
                bf128_add(bf_y[iy + 2], bf128_mul(bf_x_hat[2], bf_two));
            bf_y[iy + 2] =
                bf128_add(bf_y[iy + 2], bf128_mul(bf_x_hat[3], bf_three));

            // Step: 17
            bf_y[iy + 3] =
                bf128_add(bf_xk_hat[3], bf128_mul(bf_x_hat[0], bf_three));
            bf_y[iy + 3] = bf128_add(bf_y[iy + 3], bf_x_hat[1]);
            bf_y[iy + 3] = bf128_add(bf_y[iy + 3], bf_x_hat[2]);
            bf_y[iy + 3] =
                bf128_add(bf_y[iy + 3], bf128_mul(bf_x_hat[3], bf_two));
        }
    }
}

static void aes_enc_forward_128(const bf128_t* bf_x, const bf128_t* bf_xk,
                                const bf128_t* bf_in, const uint8_t* in,
                                uint8_t Mtag, uint8_t Mkey,
                                const uint8_t* delta, bf128_t* bf_y,
                                bool in_flag) {
    const bf128_t bf_delta = delta ? bf128_load(delta) : bf128_zero();
    const bf128_t bf_factor =
        bf128_add(bf128_mul_bit(bf_delta, Mkey), bf128_from_bit(1 ^ Mkey));

    // Step: 2..4
    for (unsigned int i = 0; i < 16; i++) {
        if (in_flag) {
            bf_y[i] = bf128_add(bf128_byte_combine(bf_in + (8 * i)),
                                bf128_byte_combine(bf_xk + (8 * i)));
        } else {
            bf128_t bf_xin[8];
            for (unsigned int j = 0; j < 8; j++) {
                bf_xin[j] =
                    bf128_mul_bit(bf_factor, (1 ^ Mtag) & get_bit(in[i], j));
            }
            // Step: 5
            bf_y[i] = bf128_add(bf128_byte_combine(bf_xin),
                                bf128_byte_combine(bf_xk + (8 * i)));
        }
    }

    const bf128_t bf_two = bf128_byte_combine_bits(2);
    const bf128_t bf_three = bf128_byte_combine_bits(3);

    for (unsigned int j = 1; j < 10; j++) {
        for (unsigned int c = 0; c <= 3; c++) {
            const unsigned int ix = 128 * (j - 1) + 32 * c;
            const unsigned int ik = 128 * j + 32 * c;
            const unsigned int iy = 16 * j + 4 * c;

            bf128_t bf_x_hat[4];
            bf128_t bf_xk_hat[4];
            for (unsigned int r = 0; r <= 3; r++) {
                // Step: 12..13
                bf_x_hat[r] = bf128_byte_combine(bf_x + (ix + 8 * r));
                bf_xk_hat[r] = bf128_byte_combine(bf_xk + (ik + 8 * r));
            }

            // Step : 14
            bf_y[iy + 0] =
                bf128_add(bf_xk_hat[0], bf128_mul(bf_x_hat[0], bf_two));
            bf_y[iy + 0] =
                bf128_add(bf_y[iy + 0], bf128_mul(bf_x_hat[1], bf_three));
            bf_y[iy + 0] = bf128_add(bf_y[iy + 0], bf_x_hat[2]);
            bf_y[iy + 0] = bf128_add(bf_y[iy + 0], bf_x_hat[3]);

            // Step: 15
            bf_y[iy + 1] = bf128_add(bf_xk_hat[1], bf_x_hat[0]);
            bf_y[iy + 1] =
                bf128_add(bf_y[iy + 1], bf128_mul(bf_x_hat[1], bf_two));
            bf_y[iy + 1] =
                bf128_add(bf_y[iy + 1], bf128_mul(bf_x_hat[2], bf_three));
            bf_y[iy + 1] = bf128_add(bf_y[iy + 1], bf_x_hat[3]);

            // Step: 16
            bf_y[iy + 2] = bf128_add(bf_xk_hat[2], bf_x_hat[0]);
            bf_y[iy + 2] = bf128_add(bf_y[iy + 2], bf_x_hat[1]);
            bf_y[iy + 2] =
                bf128_add(bf_y[iy + 2], bf128_mul(bf_x_hat[2], bf_two));
            bf_y[iy + 2] =
                bf128_add(bf_y[iy + 2], bf128_mul(bf_x_hat[3], bf_three));

            // Step: 17
            bf_y[iy + 3] =
                bf128_add(bf_xk_hat[3], bf128_mul(bf_x_hat[0], bf_three));
            bf_y[iy + 3] = bf128_add(bf_y[iy + 3], bf_x_hat[1]);
            bf_y[iy + 3] = bf128_add(bf_y[iy + 3], bf_x_hat[2]);
            bf_y[iy + 3] =
                bf128_add(bf_y[iy + 3], bf128_mul(bf_x_hat[3], bf_two));
        }
    }
}

static void aes_enc_backward_128_1(const uint8_t* x, const uint8_t* xk,
                                   const uint8_t* out, bf128_t* y_out) {
    // called only with Mtag == Mkey == 0

    uint8_t xtilde;
    // Step:2..4
    for (unsigned int j = 0; j < 10; j++) {
        for (unsigned int c = 0; c <= 3; c++) {
            for (unsigned int r = 0; r <= 3; r++) {
                // Step: 5..6
                unsigned int ird =
                    (128 * j) + (32 * ((c - r + 4) % 4)) + (8 * r);
                if (j < (10 - 1)) {
                    // Step: 7
                    xtilde = x[ird / 8];
                } else {
                    const uint8_t xout = out[(ird - 128 * (10 - 1)) / 8];
                    xtilde = xout ^ xk[(128 + ird) / 8];
                }

                const uint8_t ytilde = rotr8(xtilde, 7) ^ rotr8(xtilde, 5) ^
                                       rotr8(xtilde, 2) ^ 0x5;

                // Step: 18
                y_out[16 * j + 4 * c + r] = bf128_byte_combine_bits(ytilde);
            }
        }
    }
}

static void aes_enc_backward_128(const bf128_t* bf_x, const bf128_t* bf_xk,
                                 const bf128_t* bf_out, const uint8_t* out,
                                 uint8_t Mtag, uint8_t Mkey,
                                 const uint8_t* delta, bf128_t* y_out,
                                 bool out_flag) {
    // Step: 1
    const bf128_t bf_delta = delta ? bf128_load(delta) : bf128_zero();
    const bf128_t factor = bf128_mul_bit(
        bf128_add(bf128_mul_bit(bf_delta, Mkey), bf128_from_bit(1 ^ Mkey)),
        1 ^ Mtag);

    // Step: 2..4
    for (unsigned int j = 0; j < 10; j++) {
        for (unsigned int c = 0; c <= 3; c++) {
            for (unsigned int r = 0; r <= 3; r++) {
                bf128_t bf_x_tilde[8];
                // Step: 5
                unsigned int ird =
                    (128 * j) + (32 * ((c - r + 4) % 4)) + (8 * r);
                // Step: 6
                if (j < (10 - 1)) {
                    // Step: 7
                    memcpy(bf_x_tilde, bf_x + ird, sizeof(bf_x_tilde));
                } else {
                    // Step: 10
                    for (unsigned int i = 0; i < 8; ++i) {
                        // Step: 11
                        if (out_flag) {
                            bf_x_tilde[i] =
                                bf128_add(bf_out[ird - 128 * (10 - 1) + i],
                                          bf_xk[128 + ird + i]);
                        } else {
                            bf128_t bf_xout = bf128_mul_bit(
                                factor,
                                get_bit(out[(ird - 128 * (10 - 1)) / 8], i));
                            // Step: 12
                            bf_x_tilde[i] =
                                bf128_add(bf_xout, bf_xk[128 + ird + i]);
                        }
                    }
                }
                // Step: 13..17
                bf128_t bf_y_tilde[8];
                for (unsigned int i = 0; i < 8; ++i) {
                    bf_y_tilde[i] =
                        bf128_add(bf128_add(bf_x_tilde[(i + 7) % 8],
                                            bf_x_tilde[(i + 5) % 8]),
                                  bf_x_tilde[(i + 2) % 8]);
                }
                bf_y_tilde[0] = bf128_add(bf_y_tilde[0], factor);
                bf_y_tilde[2] = bf128_add(bf_y_tilde[2], factor);

                // Step: 18
                y_out[16 * j + 4 * c + r] = bf128_byte_combine(bf_y_tilde);
            }
        }
    }
}
// done
void aes_enc_constraints_Mkey_0_128(const uint8_t* in, const uint8_t* out,
                                    const uint8_t* w, const bf128_t* v,
                                    const uint8_t* k, const bf128_t* vk,
                                    const bf128_t* bf_in, const bf128_t* bf_out,
                                    bf128_t* a_0_vec, bf128_t* a_1_vec,
                                    bf128_t* a_2_vec, bool in_flag,
                                    bool out_flag) {
    bf128_t s[FAEST_128F_Senc];
    bf128_t vs[FAEST_128F_Senc];
    bf128_t s_dash[FAEST_128F_Senc];
    bf128_t vs_dash[FAEST_128F_Senc];
    aes_enc_forward_128_1(w, k, in, s);
    aes_enc_forward_128(v, vk, bf_in, in, 1, 0, NULL, vs, in_flag);
    aes_enc_backward_128_1(w, k, out, s_dash);
    aes_enc_backward_128(v, vk, bf_out, out, 1, 0, NULL, vs_dash, out_flag);

    for (unsigned int j = 0; j < FAEST_128F_Senc; j++) {
        bf128_t v1 = vs[j];
        bf128_t v2 = vs_dash[j];
        bf128_t u1 = s[j];
        bf128_t u2 = s_dash[j];

        // a_0_vec[j] = v1^2*v2
        bf128_t v1_squared = bf128_mul(v1, v1);
        a_0_vec[2 * j] = bf128_mul(v1_squared, v2);
        // a_1_vec[j] = v1^2*u2
        a_1_vec[2 * j] = bf128_mul(v1_squared, u2);
        // a_2_vec[j] = u1^2*v2 + v1
        bf128_t u1_squared = bf128_mul(u1, u1);
        a_2_vec[2 * j] = bf128_add(bf128_mul(u1_squared, v2), v1);

        // a_0_vec[j] = v2^2*v1
        bf128_t v2_squared = bf128_mul(v2, v2);
        a_0_vec[2 * j + 1] = bf128_mul(v2_squared, v1);
        // a_1_vec[j] = v2^2*u1
        a_1_vec[2 * j + 1] = bf128_mul(v2_squared, u1);
        // a_2_vec[j] = u2^2*v1 + v2
        bf128_t u2_squared = bf128_mul(u2, u2);
        a_2_vec[2 * j + 1] = bf128_add(bf128_mul(u2_squared, v1), v2);
    }
}

void aes_enc_constraints_Mkey_1_128(const uint8_t* in, const uint8_t* out,
                                    const bf128_t* q, const bf128_t* qk,
                                    const uint8_t* delta, const bf128_t* bf_in,
                                    const bf128_t* bf_out, bf128_t* b_vec,
                                    bool in_flag, bool out_flag) {
    // Step: 11..12
    bf128_t qs[FAEST_128F_Senc];
    bf128_t qs_dash[FAEST_128F_Senc];
    aes_enc_forward_128(q, qk, bf_in, in, 0, 1, delta, qs, in_flag);
    aes_enc_backward_128(q, qk, bf_out, out, 0, 1, delta, qs_dash, out_flag);

    // Step: 13..14
    bf128_t delta_sqr = bf128_mul(bf128_load(delta), bf128_load(delta));
    for (unsigned int j = 0; j < FAEST_128F_Senc; j++) {
        bf128_t q1_q2 = bf128_mul(qs[j], qs_dash[j]);
        b_vec[2 * j] =
            bf128_add(bf128_mul(q1_q2, qs[j]), bf128_mul(qs[j], delta_sqr));
        b_vec[2 * j + 1] = bf128_add(bf128_mul(q1_q2, qs_dash[j]),
                                     bf128_mul(qs_dash[j], delta_sqr));
    }
}

void aes_prove_128(const uint8_t* w, const bf128_t* bf_v, const uint8_t* in,
                   const uint8_t* out, bf128_t* a_0_vec, bf128_t* a_1_vec,
                   bf128_t* a_2_vec, bool key_schedule_flag, uint8_t* k,
                   bf128_t* vk, bool in_flag, bool out_flag) {
    if (key_schedule_flag) {
        aes_key_schedule_constraints_Mkey_0_128(
            w + 16, w, bf_v + 128, bf_v, a_0_vec, a_1_vec, a_2_vec, k, vk);

        aes_enc_constraints_Mkey_0_128(
            in, out, w + FAEST_128F_Lke / 8, bf_v + FAEST_128F_Lke, k, vk, NULL,
            NULL, a_0_vec + 80, a_1_vec + 80, a_2_vec + 80, in_flag, out_flag);
    } else if (in_flag && out_flag) {
        aes_enc_constraints_Mkey_0_128(in, out, w, bf_v, k, vk, bf_v + 1280,
                                       bf_v + 9 * 128, a_0_vec, a_1_vec,
                                       a_2_vec, in_flag, out_flag);
    } else if (!in_flag && out_flag) {
        aes_enc_constraints_Mkey_0_128(in, out, w, bf_v, k, vk, NULL,
                                       bf_v + 9 * 128, a_0_vec, a_1_vec,
                                       a_2_vec, in_flag, out_flag);
    }
}

void aes_verify_128(const bf128_t* bf_q, const uint8_t* delta,
                    const uint8_t* in, const uint8_t* out, bf128_t* b_vec,
                    bool key_schedule_flag, bf128_t* qk, bool in_flag,
                    bool out_flag) {
    if (key_schedule_flag) {
        aes_key_schedule_constraints_Mkey_1_128(bf_q + 128, bf_q, delta, b_vec,
                                                qk);

        aes_enc_constraints_Mkey_1_128(in, out, bf_q + FAEST_128F_Lke, qk,
                                       delta, NULL, NULL, b_vec + 80, in_flag,
                                       out_flag);
    } else if (in_flag && out_flag) {
        aes_enc_constraints_Mkey_1_128(in, out, bf_q, qk, delta, bf_q + 1280,
                                       bf_q + 9 * 128, b_vec, in_flag,
                                       out_flag);
    } else if (!in_flag && out_flag) {
        aes_enc_constraints_Mkey_1_128(in, out, bf_q, qk, delta, NULL,
                                       bf_q + 9 * 128, b_vec, in_flag,
                                       out_flag);
    }
}

static void rijnd_enc_forward_256_1(const uint8_t* x, const uint8_t* xk,
                                    const uint8_t* in, bf128_t* bf_y) {
    // Step: 2
    for (unsigned int i = 0; i < 4 * FAEST_EM_256F_Nwd; i++) {
        const uint8_t xin = in[i];
        bf_y[i] = bf128_add(bf128_byte_combine_bits(xin),
                            bf128_byte_combine_bits(xk[i]));
    }

    const bf128_t bf_two = bf128_byte_combine_bits(2);
    const bf128_t bf_three = bf128_byte_combine_bits(3);

    for (unsigned int j = 1; j < FAEST_EM_256F_R; j++) {
        for (unsigned int c = 0; c < FAEST_EM_256F_Nwd; c++) {
            const unsigned int ix = 256 * (j - 1) + 32 * c;
            const unsigned int ik = 256 * j + 32 * c;
            const unsigned int iy = 32 * j + 4 * c;

            bf128_t bf_x_hat[4];
            bf128_t bf_xk_hat[4];
            for (unsigned int r = 0; r <= 3; r++) {
                // Step: 12..13
                bf_x_hat[r] = bf128_byte_combine_bits(x[(ix + 8 * r) / 8]);
                bf_xk_hat[r] = bf128_byte_combine_bits(xk[(ik + 8 * r) / 8]);
            }

            bf_y[iy + 0] = bf128_add(bf128_mul(bf_x_hat[0], bf_two),
                                     bf128_mul(bf_x_hat[1], bf_three));
            bf_y[iy + 0] = bf128_add(bf_y[iy + 0], bf_x_hat[2]);
            bf_y[iy + 0] = bf128_add(bf_y[iy + 0], bf_x_hat[3]);
            bf_y[iy + 0] = bf128_add(bf_y[iy + 0], bf_xk_hat[0]);

            bf_y[iy + 1] =
                bf128_add(bf_x_hat[0], bf128_mul(bf_x_hat[1], bf_two));
            bf_y[iy + 1] =
                bf128_add(bf_y[iy + 1], bf128_mul(bf_x_hat[2], bf_three));
            bf_y[iy + 1] = bf128_add(bf_y[iy + 1], bf_x_hat[3]);
            bf_y[iy + 1] = bf128_add(bf_y[iy + 1], bf_xk_hat[1]);

            bf_y[iy + 2] = bf128_add(bf_x_hat[0], bf_x_hat[1]);
            bf_y[iy + 2] =
                bf128_add(bf_y[iy + 2], bf128_mul(bf_x_hat[2], bf_two));
            bf_y[iy + 2] =
                bf128_add(bf_y[iy + 2], bf128_mul(bf_x_hat[3], bf_three));
            bf_y[iy + 2] = bf128_add(bf_y[iy + 2], bf_xk_hat[2]);

            bf_y[iy + 3] =
                bf128_add(bf128_mul(bf_x_hat[0], bf_three), bf_x_hat[1]);
            bf_y[iy + 3] = bf128_add(bf_y[iy + 3], bf_x_hat[2]);
            bf_y[iy + 3] =
                bf128_add(bf_y[iy + 3], bf128_mul(bf_x_hat[3], bf_two));
            bf_y[iy + 3] = bf128_add(bf_y[iy + 3], bf_xk_hat[3]);
        }
    }
}

static void rijnd_enc_forward_256(const bf128_t* bf_x, const bf128_t* bf_xk,
                                  const bf128_t* bf_in, const uint8_t* in,
                                  uint8_t Mtag, uint8_t Mkey,
                                  const uint8_t* delta, bf128_t* bf_y,
                                  bool in_flag) {
    const bf128_t bf_delta = delta ? bf128_load(delta) : bf128_zero();
    const bf128_t bf_factor =
        bf128_add(bf128_mul_bit(bf_delta, Mkey), bf128_from_bit(1 ^ Mkey));

    for (unsigned int i = 0; i < 32; i++) {
        if (in_flag) {
            bf_y[i] = bf128_add(bf128_byte_combine(bf_in + (8 * i)),
                                bf128_byte_combine(bf_xk + (8 * i)));
        } else {
            bf128_t bf_xin[8];
            for (unsigned int j = 0; j < 8; j++) {
                bf_xin[j] =
                    bf128_mul_bit(bf_factor, (1 ^ Mtag) & get_bit(in[i], j));
            }
            // Step: 5
            bf_y[i] = bf128_add(bf128_byte_combine(bf_xin),
                                bf128_byte_combine(bf_xk + (8 * i)));
        }
    }
    const bf128_t bf_two = bf128_byte_combine_bits(2);
    const bf128_t bf_three = bf128_byte_combine_bits(3);

    for (unsigned int j = 1; j < 14; j++) {
        for (unsigned int c = 0; c < 8; c++) {
            const unsigned int ix = 256 * (j - 1) + 32 * c;
            const unsigned int ik = 256 * j + 32 * c;
            const unsigned int iy = 32 * j + 4 * c;

            bf128_t bf_x_hat[4];
            bf128_t bf_xk_hat[4];
            for (unsigned int r = 0; r <= 3; r++) {
                // Step: 12..13
                bf_x_hat[r] = bf128_byte_combine(bf_x + (ix + 8 * r));

                bf_xk_hat[r] = bf128_byte_combine(bf_xk + (ik + 8 * r));
            }

            bf_y[iy + 0] = bf128_add(bf128_mul(bf_x_hat[0], bf_two),
                                     bf128_mul(bf_x_hat[1], bf_three));
            bf_y[iy + 0] = bf128_add(bf_y[iy + 0], bf_x_hat[2]);
            bf_y[iy + 0] = bf128_add(bf_y[iy + 0], bf_x_hat[3]);
            bf_y[iy + 0] = bf128_add(bf_y[iy + 0], bf_xk_hat[0]);

            bf_y[iy + 1] =
                bf128_add(bf_x_hat[0], bf128_mul(bf_x_hat[1], bf_two));
            bf_y[iy + 1] =
                bf128_add(bf_y[iy + 1], bf128_mul(bf_x_hat[2], bf_three));
            bf_y[iy + 1] = bf128_add(bf_y[iy + 1], bf_x_hat[3]);
            bf_y[iy + 1] = bf128_add(bf_y[iy + 1], bf_xk_hat[1]);

            bf_y[iy + 2] = bf128_add(bf_x_hat[0], bf_x_hat[1]);
            bf_y[iy + 2] =
                bf128_add(bf_y[iy + 2], bf128_mul(bf_x_hat[2], bf_two));
            bf_y[iy + 2] =
                bf128_add(bf_y[iy + 2], bf128_mul(bf_x_hat[3], bf_three));
            bf_y[iy + 2] = bf128_add(bf_y[iy + 2], bf_xk_hat[2]);

            bf_y[iy + 3] =
                bf128_add(bf128_mul(bf_x_hat[0], bf_three), bf_x_hat[1]);
            bf_y[iy + 3] = bf128_add(bf_y[iy + 3], bf_x_hat[2]);
            bf_y[iy + 3] =
                bf128_add(bf_y[iy + 3], bf128_mul(bf_x_hat[3], bf_two));
            bf_y[iy + 3] = bf128_add(bf_y[iy + 3], bf_xk_hat[3]);
        }
    }
}

static void rijnd_enc_backward_256_1(const uint8_t* x, const uint8_t* xk,
                                     const uint8_t* out, bf128_t* y_out) {
    // only called with Mtag == Mkey == 0

    for (unsigned int j = 0; j < 14; j++) {
        for (unsigned int c = 0; c < 8; c++) {
            for (unsigned int r = 0; r <= 3; r++) {
                unsigned int icol = (c - r + 8) % 8;
                if (r >= 2) {
                    icol = (icol - 1 + 8) % 8;
                }
                unsigned int ird = 256 * j + 32 * icol + 8 * r;
                uint8_t x_tilde = 0;
                if (j < (14 - 1)) {
                    x_tilde = x[ird / 8];  // need to test
                } else {
                    x_tilde =
                        out[(ird - 256 * (14 - 1)) / 8] ^ xk[(ird + 256) / 8];
                }

                const uint8_t y_tilde = rotr8(x_tilde, 7) ^ rotr8(x_tilde, 5) ^
                                        rotr8(x_tilde, 2) ^ 0x5;

                // Step: 18
                y_out[32 * j + 4 * c + r] = bf128_byte_combine_bits(y_tilde);
            }
        }
    }
}

static void rijnd_enc_backward_256(const bf128_t* bf_x, const bf128_t* bf_xk,
                                   const bf128_t* bf_out, const uint8_t* out,
                                   uint8_t Mtag, uint8_t Mkey,
                                   const uint8_t* delta, bf128_t* y_out,
                                   bool out_flag) {
    // Step: 1
    const bf128_t bf_delta = delta ? bf128_load(delta) : bf128_zero();
    const bf128_t factor = bf128_mul_bit(
        bf128_add(bf128_mul_bit(bf_delta, Mkey), bf128_from_bit(1 ^ Mkey)),
        1 ^ Mtag);

    for (unsigned int j = 0; j < 14; j++) {
        for (unsigned int c = 0; c < 8; c++) {
            for (unsigned int r = 0; r <= 3; r++) {
                bf128_t bf_x_tilde[8];
                unsigned int icol = (c - r + 8) % 8;
                if (r >= 2) {
                    icol = (icol - 1 + 8) % 8;
                }
                const unsigned int ird = 256 * j + 32 * icol + 8 * r;

                if (j < (14 - 1)) {
                    memcpy(bf_x_tilde, bf_x + ird, sizeof(bf_x_tilde));
                } else {
                    for (unsigned int i = 0; i < 8; ++i) {
                        if (out_flag) {
                            bf_x_tilde[i] =
                                bf128_add(bf_out[ird - 256 * (14 - 1) + i],
                                          bf_xk[256 + ird + i]);
                        } else {
                            bf128_t bf_xout = bf128_mul_bit(
                                factor,
                                get_bit(out[(ird - 256 * (14 - 1)) / 8], i));
                            // Step: 12
                            bf_x_tilde[i] =
                                bf128_add(bf_xout, bf_xk[256 + ird + i]);
                        }
                    }
                }

                bf128_t bf_y_tilde[8];
                for (unsigned int i = 0; i < 8; ++i) {
                    bf_y_tilde[i] =
                        bf128_add(bf128_add(bf_x_tilde[(i + 7) % 8],
                                            bf_x_tilde[(i + 5) % 8]),
                                  bf_x_tilde[(i + 2) % 8]);
                }
                bf_y_tilde[0] = bf128_add(bf_y_tilde[0], factor);
                bf_y_tilde[2] = bf128_add(bf_y_tilde[2], factor);

                // Step: 18
                y_out[32 * j + 4 * c + r] = bf128_byte_combine(bf_y_tilde);
            }
        }
    }
}

// done
void rijnd_enc_constraints_Mkey_0_256(const uint8_t* in, const uint8_t* out,
                                      const uint8_t* w, const bf128_t* v,
                                      const uint8_t* k, const bf128_t* vk,
                                      const bf128_t* bf_in,
                                      const bf128_t* bf_out, bf128_t* a_0_vec,
                                      bf128_t* a_1_vec, bf128_t* a_2_vec,
                                      bool in_flag, bool out_flag) {
    bf128_t s[448];
    bf128_t vs[448];
    bf128_t s_dash[448];
    bf128_t vs_dash[448];
    rijnd_enc_forward_256_1(w, k, in, s);
    rijnd_enc_forward_256(v, vk, bf_in, in, 1, 0, NULL, vs, in_flag);
    rijnd_enc_backward_256_1(w, k, out, s_dash);
    rijnd_enc_backward_256(v, vk, bf_out, out, 1, 0, NULL, vs_dash, out_flag);

    for (unsigned int j = 0; j < 448; j++) {
        bf128_t v1 = vs[j];
        bf128_t v2 = vs_dash[j];
        bf128_t u1 = s[j];
        bf128_t u2 = s_dash[j];

        // a_0_vec[j] = v1^2*v2
        bf128_t v1_squared = bf128_mul(v1, v1);
        a_0_vec[2 * j] = bf128_mul(v1_squared, v2);
        // a_1_vec[j] = v1^2*u2
        a_1_vec[2 * j] = bf128_mul(v1_squared, u2);
        // a_2_vec[j] = u1^2*v2 + v1
        bf128_t u1_squared = bf128_mul(u1, u1);
        a_2_vec[2 * j] = bf128_add(bf128_mul(u1_squared, v2), v1);

        // a_0_vec[j] = v2^2*v1
        bf128_t v2_squared = bf128_mul(v2, v2);
        a_0_vec[2 * j + 1] = bf128_mul(v2_squared, v1);
        // a_1_vec[j] = v2^2*u1
        a_1_vec[2 * j + 1] = bf128_mul(v2_squared, u1);
        // a_2_vec[j] = u2^2*v1 + v2
        bf128_t u2_squared = bf128_mul(u2, u2);
        a_2_vec[2 * j + 1] = bf128_add(bf128_mul(u2_squared, v1), v2);
    }
}

void rijnd_enc_constraints_Mkey_1_256(const uint8_t* in, const uint8_t* out,
                                      const bf128_t* q, const bf128_t* qk,
                                      const uint8_t* delta,
                                      const bf128_t* bf_in,
                                      const bf128_t* bf_out, bf128_t* b_vec,
                                      bool in_flag, bool out_flag) {
    // Step: 11..12
    bf128_t qs[448];
    bf128_t qs_dash[448];
    rijnd_enc_forward_256(q, qk, bf_in, in, 0, 1, delta, qs, in_flag);
    rijnd_enc_backward_256(q, qk, bf_out, out, 0, 1, delta, qs_dash, out_flag);

    // Step: 13..14
    bf128_t delta_sqr = bf128_mul(bf128_load(delta), bf128_load(delta));
    for (unsigned int j = 0; j < 448; j++) {
        bf128_t q1_q2 = bf128_mul(qs[j], qs_dash[j]);
        b_vec[2 * j] =
            bf128_add(bf128_mul(q1_q2, qs[j]), bf128_mul(qs[j], delta_sqr));
        b_vec[2 * j + 1] = bf128_add(bf128_mul(q1_q2, qs_dash[j]),
                                     bf128_mul(qs_dash[j], delta_sqr));
    }
}

void rijnd_prove_256(const uint8_t* w, const bf128_t* bf_v, const uint8_t* in,
                     const uint8_t* out, bf128_t* a_0_vec, bf128_t* a_1_vec,
                     bf128_t* a_2_vec, bool key_schedule_flag, uint8_t* k,
                     bf128_t* vk, bool in_flag, bool out_flag) {
    if (key_schedule_flag) {
        rijnd_key_schedule_constraints_Mkey_0_256(
            w + 32, w, bf_v + 256, bf_v, a_0_vec, a_1_vec, a_2_vec, k, vk);

        rijnd_enc_constraints_Mkey_0_256(in, out, w + (256 + 14 * 64) / 8,
                                         bf_v + 256 + 14 * 64, k, vk, NULL,
                                         NULL, a_0_vec + 224, a_1_vec + 224,
                                         a_2_vec + 224, in_flag, out_flag);
    } else if (in_flag && out_flag) {
        rijnd_enc_constraints_Mkey_0_256(
            in, out, w, bf_v, k, vk, bf_v + 256 * 14, bf_v + 256 * 13, a_0_vec,
            a_1_vec, a_2_vec, in_flag, out_flag);
    } else if (!in_flag && out_flag) {
        rijnd_enc_constraints_Mkey_0_256(in, out, w, bf_v, k, vk, NULL,
                                         bf_v + 256 * 13, a_0_vec, a_1_vec,
                                         a_2_vec, in_flag, out_flag);
    }
}

void rijnd_verify_256(const bf128_t* bf_q, const uint8_t* delta,
                      const uint8_t* in, const uint8_t* out, bf128_t* b_vec,
                      bool key_schedule_flag, bf128_t* qk, bool in_flag,
                      bool out_flag) {
    if (key_schedule_flag) {
        rijnd_key_schedule_constraints_Mkey_1_256(bf_q + 256, bf_q, delta,
                                                  b_vec, qk);

        rijnd_enc_constraints_Mkey_1_256(in, out, bf_q + (256 + 14 * 64), qk,
                                         delta, NULL, NULL, b_vec + 224,
                                         in_flag, out_flag);
    } else if (in_flag && out_flag) {
        rijnd_enc_constraints_Mkey_1_256(in, out, bf_q, qk, delta,
                                         bf_q + 256 * 14, bf_q + 256 * 13,
                                         b_vec, in_flag, out_flag);
    } else if (!in_flag && out_flag) {
        rijnd_enc_constraints_Mkey_1_256(in, out, bf_q, qk, delta, NULL,
                                         bf_q + 256 * 13, b_vec, in_flag,
                                         out_flag);
    }
}