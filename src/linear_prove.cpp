#include "linear_prove.h"

void get_constrain_prover_linear(
    const std::vector<std::vector<uint8_t>>& matrix, field::GF2_128* v,
    field::GF2_128* A, int m) {
    for (int index = 0; index < matrix.size(); index++) {
        for (int i = 0; i < m; i++) {
            if (GET_BIT(matrix[index][i / 8], i % 8)) A[index] += v[i];
        }
    }
}

void get_constrain_verifier_linear(
    const std::vector<std::vector<uint8_t>>& matrix,
    const std::vector<uint8_t>& y, field::GF2_128* q, field::GF2_128 delta,
    field::GF2_128* B,int m) {
        for(int index=0;index<matrix.size();index++){
            for(int i=0;i<m;i++){
                if(GET_BIT(matrix[index][i/8],i%8)) B[index] += q[i];
            }
            if(y[index] == 1) B[index] += delta;
        }
    }

bool verify_constrain_linear(const field::GF2_128* A, const field::GF2_128* B, int k) {
    for(int i=0;i<k;i++){
        if(A[i]!=B[i]) return false;
    }
    return true;
}