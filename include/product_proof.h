#include <algorithm>
#include <iostream>

#include "field.h"
#include "rain.h"
#include "utils.h"

void gen_witness_dot_product(uint8_t* witness, uint8_t* vf, uint8_t* y);

void get_constrain_prover_dot_product(const field::GF2_128* v,
                                      const uint8_t* witness,
                                      field::GF2_128* A_0, field::GF2_128* A_1,
                                      int mutigate_num);

void get_constrain_verifier_dot_product(const field::GF2_128* q,
                                        field::GF2_128 delta, field::GF2_128* B,
                                        int mutigate_num);

// void gen_witness_kronecker_product(uint8_t* witness, uint8_t* b);

// void get_constrain_prover_kronecker_product(const field::GF2_128* v,
//                                             const uint8_t* witness,
//                                             field::GF2_128* A_0,
//                                             field::GF2_128* A_1,
//                                             int mutigate_num);

// void get_constrain_verifier_kronecker_product(const field::GF2_128* q,
//                                             field::GF2_128 delta,
//                                             field::GF2_128* B,
//                                             int mutigate_num);