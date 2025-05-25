#include <algorithm>
#include <cmath>
#include <iostream>

#include "field.h"
#include "utils.h"

void gen_witness(uint8_t* inputs, uint8_t* value, int m);

void FA(uint8_t* values, int base, int input_length, int index);

void HA(uint8_t* values, int base, int input_length, int index);

void get_mutigate_u(const uint8_t* witness, uint8_t* bf_y, int m);

template <typename GF2>
void get_mutigate_v(const GF2* v, GF2* bf_y, int m);

template <typename GF2>
void get_mutigate_q(const GF2* q, GF2* bf_y, int m);

template <typename GF2>
void get_constrain_prover(const GF2* v, const uint8_t* witness, GF2* A_0,
                          GF2* A_1, int mutigate_num, int m);

template <typename GF2>
void get_constrain_verifier(const GF2* q, GF2 delta, GF2* B, int mutigate_num,
                            int m);

template <typename GF2>
void convert_vec_to_field(uint8_t** vec, GF2* field_vec, const unsigned int ell,
                          const unsigned int lambda);

template <typename GF2>
void get_hamming_weight(uint8_t* witness, GF2* v, uint8_t* hamming_weight_vec,
                        GF2* hamming_weight_v, int m);

template <typename GF2>
void get_hamming_weight_q(const GF2* q, GF2* hamming_weight_q, int m);

template <typename GF2>
bool verify_hamming_weight(int hamming_weight,
                           const uint8_t* hamming_weight_vec,
                           const GF2* hamming_weight_v, const GF2* q, GF2 delta,
                           int m);