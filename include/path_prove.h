#include <iostream>

#include "field.h"
#include "rain.h"
#include "utils.h"

void rain_enc_forward_128_1(const uint8_t* witness,
                            const std::vector<uint8_t>& in,
                            field::GF2_128* bf_y);

void rain_enc_forward_128_prover(field::GF2_128* v, field::GF2_128* v_vec,
                                 field::GF2_128* bf_y);

void rain_enc_forward_128_verifier(field::GF2_128* q, field::GF2_128* q_vec,
                                   field::GF2_128 delta,
                                   const std::vector<uint8_t>& in,
                                   field::GF2_128* bf_y);

void rain_enc_backword_128_1(const uint8_t* witness, field::GF2_128* bf_y);

void rain_enc_backword_128_prover(field::GF2_128* v, field::GF2_128* bf_y);

void rain_enc_backword_128_verifier(field::GF2_128* q, field::GF2_128* bf_y);

void rain_enc_constrain_128_prover(field::GF2_128* v, field::GF2_128* v_vec,
                                   const uint8_t* witness,
                                   const std::vector<uint8_t>& in,
                                   field::GF2_128* A_0, field::GF2_128* A_1);

void rain_enc_constrain_128_verifier(field::GF2_128* q, field::GF2_128* q_vec,
                                     field::GF2_128 delta,
                                     const std::vector<uint8_t>& in,
                                     field::GF2_128* B);

void hash_forward_128_1(const uint8_t* witness, field::GF2_128* bf_y,
                        field::GF2_128* muti_gate_input);

void hash_forward_128_prover(field::GF2_128* v, field::GF2_128* v_vec,
                             field::GF2_128* bf_y,
                             field::GF2_128* muti_gate_input);

void hash_forward_128_verifier(field::GF2_128* q, field::GF2_128* q_vec,
                               field::GF2_128 delta, field::GF2_128* bf_y,
                               field::GF2_128* muti_gate_input);

void hash_backword_128_1(const uint8_t* witness, field::GF2_128* bf_y,
                         field::GF2_128* muti_gate_output);

void hash_backword_128_prover(field::GF2_128* v, field::GF2_128* bf_y,
                              field::GF2_128* muti_gate_output);

void hash_backword_128_verifier(field::GF2_128* q, field::GF2_128* bf_y,
                                field::GF2_128* muti_gate_output);

void hash_constrain_128_prover(field::GF2_128* v, field::GF2_128* v_vec,
                               const uint8_t* witness, field::GF2_128* A_0,
                               field::GF2_128* A_1);

void hash_constrain_128_verifier(field::GF2_128* q, field::GF2_128* q_vec,
                                 field::GF2_128 delta, field::GF2_128* B);

void path_prove(const uint8_t* witness, field::GF2_128* v,
                field::GF2_128* v_vec, const std::vector<uint8_t>& in,
                field::GF2_128* A_0, field::GF2_128* A_1,
                unsigned int hash_times);

void path_verify(field::GF2_128* q, field::GF2_128* q_vec, field::GF2_128 delta,
                 const std::vector<uint8_t>& in, field::GF2_128* B,
                 unsigned int hash_times);

void convert_vec_to_field(uint8_t** vec, field::GF2_128* field_vec,
                          const unsigned int ell, const unsigned int lambda);

void gen_combined_field_vec(field::GF2_128* field_vec,
                            field::GF2_128* combined_field_vec,
                            const unsigned int ell);