#include "hamming_weight_prove.h"

void gen_witness(uint8_t* inputs, uint8_t* value, int m) {
    int input_length = m;
    int base = 0;
    std::copy(inputs, inputs + (m + 7) / 8, value);
    while (input_length > 1) {
        for (int i = 0; i < (input_length - 1) / 2; i++) {
            // FA
            FA(value, base, input_length, i);
        }
        if ((input_length + 1) % 2) {
            // HA
            HA(value, base, input_length, input_length / 2 - 1);
        }
        base += input_length + input_length / 2;
        input_length /= 2;
    }
}

void FA(uint8_t* values, int base, int input_length, int index) {
    uint8_t x1, x2, x3;
    if (index == 0) {
        x1 = GET_BIT(values[base / 8], base % 8);
    } else {
        x1 = GET_BIT(values[(base + input_length + index - 1) / 8],
                     (base + input_length + index - 1) % 8);
    }

    x2 =
        GET_BIT(values[(base + 2 * index + 1) / 8], (base + 2 * index + 1) % 8);
    x3 =
        GET_BIT(values[(base + 2 * index + 2) / 8], (base + 2 * index + 2) % 8);

    setbit(&values[(base + input_length + index) / 8],
           (base + input_length + index) % 8, (x1 + x2 + x3) % 2);

    setbit(&values[(base + input_length + input_length / 2 + index) / 8],
           (base + input_length + input_length / 2 + index) % 8,
           (((x1 + x2) % 2) * ((x1 + x3) % 2) + x1) % 2);
}

void HA(uint8_t* values, int base, int input_length, int index) {
    uint8_t x1, x2;
    if (index == 0) {
        x1 = GET_BIT(values[base / 8], base % 8);
    } else {
        x1 = GET_BIT(values[(base + input_length + index - 1) / 8],
                     (base + input_length + index - 1) % 8);
    }

    x2 = GET_BIT(values[(base + input_length - 1) / 8],
                 (base + input_length - 1) % 8);

    setbit(&values[(base + input_length + input_length / 2 - 1) / 8],
           (base + input_length + input_length / 2 - 1) % 8, (x1 + x2) % 2);
    setbit(&values[(base + input_length + input_length / 2 + input_length / 2 -
                    1) /
                   8],
           (base + input_length + input_length / 2 + input_length / 2 - 1) % 8,
           x1 * x2);
}

void get_mutigate_u(const uint8_t* witness, uint8_t* bf_y, int m) {
    int input_length = m;
    int base = 0;
    int index = 0;
    uint8_t x1, x2, x3, x4;
    while (input_length > 1) {
        for (int i = 0; i < (input_length - 1) / 2; i++) {
            // FA
            if (i == 0) {
                x1 = GET_BIT(witness[base / 8], base % 8);
            } else {
                x1 = GET_BIT(witness[(base + input_length + i - 1) / 8],
                             (base + input_length + i - 1) % 8);
            }
            x2 = GET_BIT(witness[(base + 2 * i + 1) / 8],
                         (base + 2 * i + 1) % 8);
            x3 = GET_BIT(witness[(base + 2 * i + 2) / 8],
                         (base + 2 * i + 2) % 8);

            x4 = GET_BIT(
                witness[(base + input_length + input_length / 2 + i) / 8],
                (base + input_length + input_length / 2 + i) % 8);

            setbit(&bf_y[index], 0, (x1 + x2) % 2);
            setbit(&bf_y[index], 1, (x1 + x3) % 2);
            setbit(&bf_y[index], 2, (x1 + x4) % 2);
            index++;
        }
        if ((input_length + 1) % 2) {
            // HA
            if (input_length == 2) {
                x1 = GET_BIT(witness[base / 8], base % 8);
            } else {
                x1 = GET_BIT(
                    witness[(base + input_length + input_length / 2 - 2) / 8],
                    (base + input_length + input_length / 2 - 2) % 8);
            }
            x2 = GET_BIT(witness[(base + input_length - 1) / 8],
                         (base + input_length - 1) % 8);
            setbit(&bf_y[index], 0, x1);
            setbit(&bf_y[index], 1, x2);
            setbit(&bf_y[index], 2, x1 * x2);
            index++;
        }
        base += input_length + input_length / 2;
        input_length /= 2;
    }
}


template <typename GF2>
void get_mutigate_v(const GF2* v, GF2* bf_y, int m) {
    int input_length = m;
    int base = 0;
    int index = 0;
    GF2 v1, v2, v3, v4;
    while (input_length > 1) {
        for (int i = 0; i < (input_length - 1) / 2; i++) {
            // FA
            if (i == 0) {
                v1 = v[base];
            } else {
                v1 = v[base + input_length + i - 1];
            }
            v2 = v[base + 2 * i + 1];
            v3 = v[base + 2 * i + 2];
            v4 = v[base + input_length + input_length / 2 + i];
            bf_y[index++] = v1 + v2;
            bf_y[index++] = v1 + v3;
            bf_y[index++] = v1 + v4;
        }
        if ((input_length + 1) % 2) {
            // HA
            if (input_length == 2) {
                v1 = v[base];
            } else {
                v1 = v[base + input_length + input_length / 2 - 2];
            }
            v2 = v[base + input_length - 1];
            v4 = v[base + input_length + input_length / 2 + input_length / 2 -
                   1];
            bf_y[index++] = v1;
            bf_y[index++] = v2;
            bf_y[index++] = v4;
        }
        base += input_length + input_length / 2;
        input_length /= 2;
    }
}


template <typename GF2>
void get_mutigate_q(const GF2* q, GF2* bf_y, int m) {
    int input_length = m;
    int base = 0;
    int index = 0;
    GF2 q1, q2, q3, q4;
    while (input_length > 1) {
        for (int i = 0; i < (input_length - 1) / 2; i++) {
            // FA
            if (i == 0) {
                q1 = q[base];
            } else {
                q1 = q[base + input_length + i - 1];
            }
            q2 = q[base + 2 * i + 1];
            q3 = q[base + 2 * i + 2];
            q4 = q[base + input_length + input_length / 2 + i];
            bf_y[index++] = q1 + q2;
            bf_y[index++] = q1 + q3;
            bf_y[index++] = q1 + q4;
        }
        if ((input_length + 1) % 2) {
            // HA
            if (input_length == 2) {
                q1 = q[base];
            } else {
                q1 = q[base + input_length + input_length / 2 - 2];
            }
            q2 = q[base + input_length - 1];
            q4 = q[base + input_length + input_length / 2 + input_length / 2 -
                   1];
            bf_y[index++] = q1;
            bf_y[index++] = q2;
            bf_y[index++] = q4;
        }
        base += input_length + input_length / 2;
        input_length /= 2;
    }
}

template <typename GF2>
void get_constrain_prover(const GF2* v, const uint8_t* witness,
                          GF2* A_0, GF2* A_1,
                          int mutigate_num, int m) {
    std::vector<uint8_t> mutigate_u(mutigate_num);
    std::vector<GF2> mutigate_v(3 * mutigate_num);
    get_mutigate_u(witness, mutigate_u.data(), m);
    get_mutigate_v(v, mutigate_v.data(), m);
    for (int i = 0; i < mutigate_num; i++) {
        A_0[i] = mutigate_v[3 * i] * mutigate_v[3 * i + 1];
        if (GET_BIT(mutigate_u[i], 0) == 1) A_1[i] += mutigate_v[3 * i + 1];
        if (GET_BIT(mutigate_u[i], 1) == 1) A_1[i] += mutigate_v[3 * i];
        A_1[i] -= mutigate_v[3 * i + 2];
    }
}

template <typename GF2>
void get_constrain_verifier(const GF2* q, GF2 delta,
                            GF2* B, int mutigate_num, int m) {
    std::vector<GF2> mutigate_q(3 * mutigate_num);
    get_mutigate_q(q, mutigate_q.data(), m);
    for (int i = 0; i < mutigate_num; i++) {
        B[i] = mutigate_q[3 * i] * mutigate_q[3 * i + 1] -
               mutigate_q[3 * i + 2] * delta;
    }
}

template <typename GF2>
void convert_vec_to_field(uint8_t** vec, GF2* field_vec,
                          const unsigned int ell, const unsigned int lambda) {
    for (unsigned int row = 0; row < ell; row++) {
        uint8_t new_row[lambda / 8] = {0};
        for (unsigned int column = 0; column < lambda; ++column) {
            ptr_set_bit(new_row, ptr_get_bit(vec[column], row), column);
        }
        field_vec[row / 8 * 8 + (row % 8)].from_bytes(new_row);
    }
}

template <typename GF2>
void get_hamming_weight(uint8_t* witness, GF2* v,
                        uint8_t* hamming_weight_vec,
                        GF2* hamming_weight_v, int m) {
    int input_length = m;
    int base = 0;
    int index = 0;
    while (input_length > 1) {
        hamming_weight_vec[index] =
            GET_BIT(witness[(base + input_length + input_length / 2 - 1) / 8],
                    (base + input_length + input_length / 2 - 1) % 8);
        hamming_weight_v[index++] =
            v[base + input_length + input_length / 2 - 1];
        base += input_length + input_length / 2;
        input_length /= 2;
    }
}

template <typename GF2>
void get_hamming_weight_q(const GF2* q, GF2* hamming_weight_q,
                          int m) {
    int input_length = m;
    int base = 0;
    int index = 0;
    while (input_length > 1) {
        hamming_weight_q[index++] =
            q[base + input_length + input_length / 2 - 1];
        base += input_length + input_length / 2;
        input_length /= 2;
    }
}

template <typename GF2>
bool verify_hamming_weight(int hamming_weight, const uint8_t* hamming_weight_vec,
                           const GF2* hamming_weight_v, const GF2* q,
                           GF2 delta, int m) {
    int vec_length = std::floor(log2(m) + 1e-10);
    std::vector<GF2> hamming_weight_q(vec_length);
    int cnt = 0;
    get_hamming_weight_q(q,hamming_weight_q.data(),m);

    for (int i = vec_length - 1; i >= 0; i--) {
        if (hamming_weight_vec[i] == 1) {
            cnt += (1 << i);
            if (hamming_weight_v[i] + delta != hamming_weight_q[i]){
                return false;
            }
         }else{
                if (hamming_weight_v[i] != hamming_weight_q[i]){
                return false;
            }
         }
    }
    return cnt == hamming_weight;
}






template void get_mutigate_v<field::GF2_128>(const field::GF2_128* v, field::GF2_128* bf_y, int m);

template void get_mutigate_q<field::GF2_128>(const field::GF2_128* q, field::GF2_128* bf_y, int m);

template void get_constrain_prover<field::GF2_128>(const field::GF2_128* v, const uint8_t* witness, field::GF2_128* A_0,
                          field::GF2_128* A_1, int mutigate_num, int m);

template void get_constrain_verifier<field::GF2_128>(const field::GF2_128* q, field::GF2_128 delta, field::GF2_128* B, int mutigate_num,
                            int m);

template void convert_vec_to_field<field::GF2_128>(uint8_t** vec, field::GF2_128* field_vec, const unsigned int ell,
                          const unsigned int lambda);

template void get_hamming_weight<field::GF2_128>(uint8_t* witness, field::GF2_128* v, uint8_t* hamming_weight_vec,
                        field::GF2_128* hamming_weight_v, int m);

template void get_hamming_weight_q<field::GF2_128>(const field::GF2_128* q, field::GF2_128* hamming_weight_q, int m);

template bool verify_hamming_weight<field::GF2_128>(int hamming_weight,
                           const uint8_t* hamming_weight_vec,
                           const field::GF2_128* hamming_weight_v, const field::GF2_128* q, field::GF2_128 delta,
                           int m);




template void get_mutigate_v<field::GF2_256>(const field::GF2_256* v, field::GF2_256* bf_y, int m);

template void get_mutigate_q<field::GF2_256>(const field::GF2_256* q, field::GF2_256* bf_y, int m);

template void get_constrain_prover<field::GF2_256>(const field::GF2_256* v, const uint8_t* witness, field::GF2_256* A_0,
                          field::GF2_256* A_1, int mutigate_num, int m);

template void get_constrain_verifier<field::GF2_256>(const field::GF2_256* q, field::GF2_256 delta, field::GF2_256* B, int mutigate_num,
                            int m);

template void convert_vec_to_field<field::GF2_256>(uint8_t** vec, field::GF2_256* field_vec, const unsigned int ell,
                          const unsigned int lambda);

template void get_hamming_weight<field::GF2_256>(uint8_t* witness, field::GF2_256* v, uint8_t* hamming_weight_vec,
                        field::GF2_256* hamming_weight_v, int m);

template void get_hamming_weight_q<field::GF2_256>(const field::GF2_256* q, field::GF2_256* hamming_weight_q, int m);

template bool verify_hamming_weight<field::GF2_256>(int hamming_weight,
                           const uint8_t* hamming_weight_vec,
                           const field::GF2_256* hamming_weight_v, const field::GF2_256* q, field::GF2_256 delta,
                           int m);