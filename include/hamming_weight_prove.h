#include <algorithm>
#include <cmath>
#include <iostream>

#include "field.h"
#include "utils.h"

void gen_witness(uint8_t* inputs, uint8_t* value, int m);

void FA(uint8_t* values, int base, int input_length, int index);

void HA(uint8_t* values, int base, int input_length, int index);

void get_mutigate_u(const uint8_t* witness, uint8_t* bf_y, int m);

void get_mutigate_v(const field::GF2_128* v, field::GF2_128* bf_y, int m);

void get_mutigate_q(const field::GF2_128* q, field::GF2_128* bf_y, int m);

void get_constrain_prover(const field::GF2_128* v, const uint8_t* witness,
                          field::GF2_128* A_0, field::GF2_128* A_1,
                          int mutigate_num, int m);

void get_constrain_verifier(const field::GF2_128* q, field::GF2_128 delta,
                            field::GF2_128* B, int mutigate_num, int m);

void convert_vec_to_field(uint8_t** vec, field::GF2_128* field_vec,
                          const unsigned int ell, const unsigned int lambda);

void get_hamming_weight(uint8_t* witness, field::GF2_128* v,
                        uint8_t* hamming_weight_vec,
                        field::GF2_128* hamming_weight_v, int m);

void get_hamming_weight_q(const field::GF2_128* q, field::GF2_128* hamming_weight_q,
                          int m);


bool verify_hamming_weight(int hamming_weight, const uint8_t* hamming_weight_vec,
                           const field::GF2_128* hamming_weight_v, const field::GF2_128* q,
                           field::GF2_128 delta, int m);