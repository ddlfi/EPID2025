#include <algorithm>
#include <cmath>
#include <iostream>

#include "field.h"
#include "utils.h"

template <typename GF2>
void get_constrain_prover_linear(
    const std::vector<std::vector<uint8_t>>& matrix, GF2* v,
    GF2* A, int m) {
    for (int index = 0; index < matrix.size(); index++) {
        for (int i = 0; i < m; i++) {
            if (GET_BIT(matrix[index][i / 8], i % 8)) A[index] += v[i];
        }
    }
}
template <typename GF2>
void get_constrain_verifier_linear(
    const std::vector<std::vector<uint8_t>>& matrix,
    const std::vector<uint8_t>& y, GF2* q, GF2 delta,
    GF2* B,int m) {
        for(int index=0;index<matrix.size();index++){
            for(int i=0;i<m;i++){
                if(GET_BIT(matrix[index][i/8],i%8)) B[index] += q[i];
            }
            if(y[index] == 1) B[index] += delta;
        }
    }
template <typename GF2>
bool verify_constrain_linear(const GF2* A, const GF2* B, int k) {
    for(int i=0;i<k;i++){
        if(A[i]!=B[i]) return false;
    }
    return true;
}