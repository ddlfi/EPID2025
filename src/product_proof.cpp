#include "product_proof.h"

void gen_witness_dot_product(uint8_t* witness, uint8_t* vf, uint8_t* y) {
    std::copy(vf, vf + 32, witness);
    std::copy(y, y + 32, witness + 32);
    for (int i = 0; i < 32; i++) {
        witness[64 + i] = vf[i] & y[i];
    }
}

void get_constrain_prover_dot_product(const field::GF2_128* v,
                                      const uint8_t* witness,
                                      field::GF2_128* A_0, field::GF2_128* A_1,
                                      int mutigate_num) {
    for (int i = 0; i < mutigate_num; i++) {
        A_0[i] = v[i] * v[i + mutigate_num];
        if (GET_BIT(witness[i / 8], i % 8) == 1) A_1[i] += v[i + mutigate_num];
        if (GET_BIT(witness[i / 8 + 32], i % 8) == 1) A_1[i] += v[i];
        A_1[i] -= v[i + 2 * mutigate_num];
    }
}

void get_constrain_verifier_dot_product(const field::GF2_128* q,
                                        field::GF2_128 delta, field::GF2_128* B,
                                        int mutigate_num) {
    for (int i = 0; i < mutigate_num; i++) {
        B[i] = q[i] * q[i + mutigate_num] - q[i + 2 * mutigate_num] * delta;
    }
}