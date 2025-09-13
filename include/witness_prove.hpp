#ifndef WITNESS_PROVE_HPP
#define WITNESS_PROVE_HPP

#include <stdint.h>
#include <vector>

#include "fields.h"

// Forward declaration
struct signature_t;

void witness_prove(const uint8_t* w, const bf128_t* bf_v, bf128_t* a_0_vec,
                   bf128_t* a_1_vec, bf128_t* a_2_vec, uint8_t* chall_2, unsigned int lambda,
                   uint8_t* t, uint8_t* r, unsigned int height,
                   unsigned int ell_hat, const std::vector<signature_t>& srl);

void witness_verify(const bf128_t* bf_q, const uint8_t* delta, bf128_t* b_vec,
                    unsigned int lambda, const uint8_t* t, const uint8_t* r,
                    unsigned int height, unsigned int ell_hat,
                    const std::vector<signature_t>& srl);

void or_prove(const uint8_t* w, const bf128_t* bf_v, bf128_t* a_0_vec,
              bf128_t* a_1_vec, bf128_t* a_2_vec, unsigned int lambda, const uint8_t* r_j,
              const uint8_t* t_j);

void or_verify(const bf128_t* bf_q, const uint8_t* delta, bf128_t* b_vec,
               unsigned int lambda, const uint8_t* r_j, const uint8_t* t_j);

void leave_prove_128(const uint8_t* w, const bf128_t* bf_v, bf128_t* a_0_vec,
                     bf128_t* a_1_vec, bf128_t* a_2_vec);

void leave_verify_128(const bf128_t* bf_q, const uint8_t* delta, bf128_t* b_vec);

void merkle_tree_prove_128(const uint8_t* w, const bf128_t* bf_v, bf128_t* a_0_vec,
                           bf128_t* a_1_vec, bf128_t* a_2_vec, unsigned int height);

void merkle_tree_verify_128(const bf128_t* bf_q, const uint8_t* delta, bf128_t* b_vec,
                            unsigned int height);

#endif  // WITNESS_PROVE_HPP