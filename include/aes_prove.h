#ifndef AES_PROVE_H
#define AES_PROVE_H

#include <stdint.h>

#include "fields.h"
#ifdef __cplusplus
extern "C" {
#endif

bf128_t* column_to_row_major_and_shrink_V_128(uint8_t** v, unsigned int ell);
bf256_t* column_to_row_major_and_shrink_V_256(uint8_t** v, unsigned int ell);

void aes_prove_128(const uint8_t* w, const bf128_t* bf_v,
                   const uint8_t* in, const uint8_t* out,
                   bf128_t* a_0_vec, bf128_t* a_1_vec, bf128_t* a_2_vec,
                   bool key_schedule_flag, uint8_t* k, bf128_t* vk,
                   bool in_flag, bool out_flag);

void aes_verify_128(const bf128_t* bf_q, const uint8_t* delta,
                    const uint8_t* in, const uint8_t* out,
                    bf128_t* b_vec, bool key_schedule_flag, bf128_t* qk,
                    bool in_flag, bool out_flag);

void rijnd_prove_256(const uint8_t* w, const bf128_t* bf_v,
                     const uint8_t* in, const uint8_t* out,
                     bf128_t* a_0_vec, bf128_t* a_1_vec, bf128_t* a_2_vec,
                     bool key_schedule_flag, uint8_t* k, bf128_t* vk,
                     bool in_flag, bool out_flag);

void rijnd_verify_256(const bf128_t* bf_q, const uint8_t* delta,
                      const uint8_t* in, const uint8_t* out,
                      bf128_t* b_vec, bool key_schedule_flag, bf128_t* qk,
                      bool in_flag, bool out_flag);

void aes_key_schedule_constraints_Mkey_0_128(
    const uint8_t* w, const uint8_t* origin_key, const bf128_t* v,
    const bf128_t* origin_key_v, bf128_t* a_0_vec, bf128_t* a_1_vec, bf128_t* a_2_vec, uint8_t* k,
    bf128_t* vk);

void aes_key_schedule_constraints_Mkey_1_128(const bf128_t* q,
                                             const bf128_t* origin_key_q,
                                             const uint8_t* delta,
                                             bf128_t* b_vec,
                                             bf128_t* qk);

void rijnd_key_schedule_constraints_Mkey_0_256(
    const uint8_t* w, const uint8_t* origin_key, const bf128_t* v,
    const bf128_t* origin_key_v, bf128_t* a_0_vec, bf128_t* a_1_vec, bf128_t* a_2_vec, uint8_t* k,
    bf128_t* vk);

void rijnd_key_schedule_constraints_Mkey_1_256(const bf128_t* q,
                                               const bf128_t* origin_key_q,
                                               const uint8_t* delta,
                                               bf128_t* b_vec,
                                               bf128_t* qk);

void aes_enc_constraints_Mkey_0_128(
    const uint8_t* in, const uint8_t* out, const uint8_t* w, const bf128_t* v,
    const uint8_t* k, const bf128_t* vk, const bf128_t* bf_in,
    const bf128_t* bf_out, bf128_t* a_0_vec, bf128_t* a_1_vec, bf128_t* a_2_vec, bool in_flag,
    bool out_flag);

void aes_enc_constraints_Mkey_1_128(
    const uint8_t* in, const uint8_t* out, const bf128_t* q, const bf128_t* qk,
    const uint8_t* delta, const bf128_t* bf_in, const bf128_t* bf_out,
    bf128_t* b_vec, bool in_flag, bool out_flag);

void rijnd_enc_constraints_Mkey_0_256(
    const uint8_t* in, const uint8_t* out, const uint8_t* w, const bf128_t* v,
    const uint8_t* k, const bf128_t* vk, const bf128_t* bf_in,
    const bf128_t* bf_out, bf128_t* a_0_vec, bf128_t* a_1_vec, bf128_t* a_2_vec, bool in_flag,
    bool out_flag);

void rijnd_enc_constraints_Mkey_1_256(
    const uint8_t* in, const uint8_t* out, const bf128_t* q, const bf128_t* qk,
    const uint8_t* delta, const bf128_t* bf_in, const bf128_t* bf_out,
    bf128_t* b_vec, bool in_flag, bool out_flag);

#ifdef __cplusplus
}
#endif

#endif  // AES_PROVE_H