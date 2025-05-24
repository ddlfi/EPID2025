#include <algorithm>
#include <cmath>
#include <iostream>

#include "field.h"
#include "utils.h"

void get_constrain_prover_linear(
    const std::vector<std::vector<uint8_t>>& matrix, field::GF2_128* v,
    field::GF2_128* A, int m);

void get_constrain_verifier_linear(
    const std::vector<std::vector<uint8_t>>& matrix,
    const std::vector<uint8_t>& y, field::GF2_128* q, field::GF2_128 delta,
    field::GF2_128* B,int m);

bool verify_constrain_linear(const field::GF2_128* A, const field::GF2_128* B, int k);