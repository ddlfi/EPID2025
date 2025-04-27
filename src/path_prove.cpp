#include "path_prove.h"

void rain_enc_forward_128_1(const uint8_t* witness,
                            const std::vector<uint8_t>& in,
                            field::GF2_128* bf_y) {
    // case: not vole,the real value
    field::GF2_128 key;
    key.from_bytes(witness);
    witness += 2 * sizeof(uint64_t);
    field::GF2_128 inverse_input_value;
    inverse_input_value.from_bytes(in.data());
    inverse_input_value += key;
    inverse_input_value += roundconst[0];
    bf_y[0] = inverse_input_value;

    for (int i = 1; i < 3; i++) {
        field::GF2_128 tmp;
        tmp.from_bytes(witness);
        witness += 2 * sizeof(uint64_t);
        tmp = tmp.multiply_with_transposed_GF2_matrix(matrix_transposed[i - 1]);
        tmp += key;
        tmp += roundconst[i];
        bf_y[i] = tmp;
    }
}

void rain_enc_forward_128_prover(field::GF2_128* v, field::GF2_128* v_vec,
                                 field::GF2_128* bf_y) {
    bf_y[0] = v[0];
    for (int i = 1; i < 3; i++) {
        bf_y[i] = gf128_vec_muti_with_transposed_GF2_matrix(
            v_vec + 128 * i, matrix_transposed[i - 1]);
        bf_y[i] += v[0];
    }
}

void rain_enc_forward_128_verifier(field::GF2_128* q, field::GF2_128* q_vec,
                                   field::GF2_128 delta,
                                   const std::vector<uint8_t>& in,
                                   field::GF2_128* bf_y) {
    field::GF2_128 in_128;
    in_128.from_bytes(in.data());
    bf_y[0] = q[0] + delta * roundconst[0] + delta * in_128;
    for (int i = 1; i < 3; i++) {
        bf_y[i] = gf128_vec_muti_with_transposed_GF2_matrix(
            q_vec + 128 * i, matrix_transposed[i - 1]);
        bf_y[i] += delta * roundconst[i] + q[0];
    }
}

void rain_enc_backword_128_1(const uint8_t* witness, field::GF2_128* bf_y) {
    field::GF2_128 key;
    key.from_bytes(witness);
    witness += 2 * sizeof(uint64_t);
    field::GF2_128 tmp;
    for (int i = 0; i < 2; i++) {
        tmp.from_bytes(witness);
        witness += 2 * sizeof(uint64_t);
        bf_y[i] = tmp;
    }

    tmp.from_bytes(witness);

    tmp += key;
    bf_y[2] = tmp;
}

void rain_enc_backword_128_prover(field::GF2_128* v, field::GF2_128* bf_y) {
    for (int i = 0; i < 2; i++) {
        bf_y[i] = v[i + 1];
    }
    bf_y[2] = v[0] + v[3];
}

void rain_enc_backword_128_verifier(field::GF2_128* q, field::GF2_128* bf_y) {
    for (int i = 0; i < 2; i++) {
        bf_y[i] = q[i + 1];
    }
    bf_y[2] = q[0] + q[3];
}

void rain_enc_constrain_128_prover(field::GF2_128* v, field::GF2_128* v_vec,
                                   const uint8_t* witness,
                                   const std::vector<uint8_t>& in,
                                   field::GF2_128* A_0, field::GF2_128* A_1) {
    std::vector<field::GF2_128> inverse_input_real_value(3);
    std::vector<field::GF2_128> inverse_output_real_value(3);
    std::vector<field::GF2_128> inverse_input_vole_value(3);
    std::vector<field::GF2_128> inverse_output_vole_value(3);
    rain_enc_forward_128_1(witness, in, inverse_input_real_value.data());
    rain_enc_forward_128_prover(v, v_vec, inverse_input_vole_value.data());
    rain_enc_backword_128_1(witness, inverse_output_real_value.data());
    rain_enc_backword_128_prover(v, inverse_output_vole_value.data());
    for (int i = 0; i < 3; i++) {
        A_0[i] = inverse_input_vole_value[i] * inverse_output_vole_value[i];
        A_1[i] = inverse_input_vole_value[i] * inverse_output_real_value[i] +
                 inverse_output_vole_value[i] * inverse_input_real_value[i];
    }
}

void rain_enc_constrain_128_verifier(field::GF2_128* q, field::GF2_128* q_vec,
                                     field::GF2_128 delta,
                                     const std::vector<uint8_t>& in,
                                     field::GF2_128* B) {
    std::vector<field::GF2_128> inverse_input_vole_value(3);
    std::vector<field::GF2_128> inverse_output_vole_value(3);
    rain_enc_forward_128_verifier(q, q_vec, delta, in,
                                  inverse_input_vole_value.data());
    rain_enc_backword_128_verifier(q, inverse_output_vole_value.data());
    for (int i = 0; i < 3; i++) {
        B[i] = inverse_input_vole_value[i] * inverse_output_vole_value[i] -
               delta * delta;
    }
}

void hash_forward_128_1(const uint8_t* witness, field::GF2_128* bf_y,
                        field::GF2_128* muti_gate_input) {
    field::GF2_128 i_0, i_1, s;
    i_0.from_bytes(witness);
    i_1.from_bytes(witness + 2 * sizeof(uint64_t));
    s.from_bytes(witness + 4 * sizeof(uint64_t));
    witness += 6 * sizeof(uint64_t);

    muti_gate_input[0] = s;
    muti_gate_input[1] = i_0 + i_1;

    field::GF2_128 s_;
    s_.from_bytes(witness);
    witness += 2 * sizeof(uint64_t);
    field::GF2_128 key, msg;
    key = s_ + i_1;
    msg = s_ + i_0;

    bf_y[0] = key + roundconst[0] + msg;
    for (int i = 1; i < NUM_SBOXES; i++) {
        field::GF2_128 tmp;
        tmp.from_bytes(witness);
        witness += 2 * sizeof(uint64_t);
        tmp = tmp.multiply_with_transposed_GF2_matrix(matrix_transposed[i - 1]);
        tmp += key + roundconst[i];
        bf_y[i] = tmp;
    }
}

void hash_forward_128_prover(field::GF2_128* v, field::GF2_128* v_vec,
                             field::GF2_128* bf_y,
                             field::GF2_128* muti_gate_input) {
    muti_gate_input[0] = v[2];
    muti_gate_input[1] = v[0] + v[1];

    field::GF2_128 v_key, v_msg;
    v_key = v[3] + v[1];
    v_msg = v[3] + v[0];
    bf_y[0] = v_key + v_msg;

    for (int i = 1; i < NUM_SBOXES; i++) {
        bf_y[i] = gf128_vec_muti_with_transposed_GF2_matrix(
            v_vec + (i + 3) * 128, matrix_transposed[i - 1]);
        bf_y[i] += v_key;
    }
}

void hash_forward_128_verifier(field::GF2_128* q, field::GF2_128* q_vec,
                               field::GF2_128 delta, field::GF2_128* bf_y,
                               field::GF2_128* muti_gate_input) {
    muti_gate_input[0] = q[2];
    muti_gate_input[1] = q[0] + q[1];

    field::GF2_128 q_key, q_msg;
    q_key = q[3] + q[1];
    q_msg = q[3] + q[0];
    bf_y[0] = q_msg + q_key + delta * roundconst[0];

    for (int i = 1; i < NUM_SBOXES; i++) {
        bf_y[i] = gf128_vec_muti_with_transposed_GF2_matrix(
            q_vec + (i + 3) * 128, matrix_transposed[i - 1]);
        bf_y[i] += q_key + delta * roundconst[i];
    }
}
// //RAIN3
// void hash_backword_128_1(const uint8_t* witness, field::GF2_128* bf_y,
//                          field::GF2_128* muti_gate_output) {
//     field::GF2_128 i_0, i_1;
//     i_0.from_bytes(witness);
//     i_1.from_bytes(witness + 2 * sizeof(uint64_t));

//     field::GF2_128 s_, inver_out_1, inver_out_2, y_out;
//     s_.from_bytes(witness + 6 * sizeof(uint64_t));
//     inver_out_1.from_bytes(witness + 8 * sizeof(uint64_t));
//     inver_out_2.from_bytes(witness + 10 * sizeof(uint64_t));
//     y_out.from_bytes(witness + 12 * sizeof(uint64_t));

//     muti_gate_output[0] = s_;
//     bf_y[0] = inver_out_1;
//     bf_y[1] = inver_out_2;
//     bf_y[2] = y_out + i_0 + i_1;
// }

// void hash_backword_128_prover(field::GF2_128* v, field::GF2_128* bf_y,
//                               field::GF2_128* muti_gate_output) {
//     muti_gate_output[0] = v[3];
//     bf_y[0] = v[4];
//     bf_y[1] = v[5];
//     bf_y[2] = v[0] + v[1] + v[6];
// }

// void hash_backword_128_verifier(field::GF2_128* q, field::GF2_128* bf_y,
//                                 field::GF2_128* muti_gate_output) {
//     muti_gate_output[0] = q[3];
//     bf_y[0] = q[4];
//     bf_y[1] = q[5];
//     bf_y[2] = q[0] + q[1] + q[6];
// }

void hash_backword_128_1(const uint8_t* witness, field::GF2_128* bf_y,
                         field::GF2_128* muti_gate_output) {
    field::GF2_128 i_0, i_1;
    i_0.from_bytes(witness);
    i_1.from_bytes(witness + 2 * sizeof(uint64_t));

    field::GF2_128 s_, inver_out_1, inver_out_2,inver_out_3, y_out;
    s_.from_bytes(witness + 6 * sizeof(uint64_t));
    inver_out_1.from_bytes(witness + 8 * sizeof(uint64_t));
    inver_out_2.from_bytes(witness + 10 * sizeof(uint64_t));
    inver_out_3.from_bytes(witness + 12 * sizeof(uint64_t));
    y_out.from_bytes(witness + 14 * sizeof(uint64_t));

    muti_gate_output[0] = s_;
    bf_y[0] = inver_out_1;
    bf_y[1] = inver_out_2;
    bf_y[2] = inver_out_3;
    bf_y[3] = y_out + i_0 + i_1;
}

void hash_backword_128_prover(field::GF2_128* v, field::GF2_128* bf_y,
                              field::GF2_128* muti_gate_output) {
    muti_gate_output[0] = v[3];
    bf_y[0] = v[4];
    bf_y[1] = v[5];
    bf_y[2] = v[6];
    bf_y[3] = v[0] + v[1] + v[7];
}

void hash_backword_128_verifier(field::GF2_128* q, field::GF2_128* bf_y,
                                field::GF2_128* muti_gate_output) {
    muti_gate_output[0] = q[3];
    bf_y[0] = q[4];
    bf_y[1] = q[5];
    bf_y[2] = q[6];
    bf_y[3] = q[0] + q[1] + q[7];
}

void hash_constrain_128_prover(field::GF2_128* v, field::GF2_128* v_vec,
                               const uint8_t* witness, field::GF2_128* A_0,
                               field::GF2_128* A_1) {
    std::vector<field::GF2_128> inverse_input_real_value(NUM_SBOXES);
    std::vector<field::GF2_128> inverse_output_real_value(NUM_SBOXES);
    std::vector<field::GF2_128> inverse_input_vole_value(NUM_SBOXES);
    std::vector<field::GF2_128> inverse_output_vole_value(NUM_SBOXES);
    std::vector<field::GF2_128> muti_gate_real_value(3);  // input:0,1 output:2
    std::vector<field::GF2_128> muti_gate_vole_value(3);

    hash_forward_128_1(witness, inverse_input_real_value.data(),
                       muti_gate_real_value.data());
    hash_forward_128_prover(v, v_vec, inverse_input_vole_value.data(),
                            muti_gate_vole_value.data());
    hash_backword_128_1(witness, inverse_output_real_value.data(),
                        muti_gate_real_value.data() + 2);
    hash_backword_128_prover(v, inverse_output_vole_value.data(),
                             muti_gate_vole_value.data() + 2);

    A_0[0] = muti_gate_vole_value[0] * muti_gate_vole_value[1];
    A_1[0] = muti_gate_vole_value[0] * muti_gate_real_value[1] +
             muti_gate_vole_value[1] * muti_gate_real_value[0] -
             muti_gate_vole_value[2];
    for (int i = 0; i < NUM_SBOXES; i++) {
        A_0[i + 1] = inverse_input_vole_value[i] * inverse_output_vole_value[i];
        A_1[i + 1] =
            inverse_input_vole_value[i] * inverse_output_real_value[i] +
            inverse_output_vole_value[i] * inverse_input_real_value[i];
    }
}

void hash_constrain_128_verifier(field::GF2_128* q, field::GF2_128* q_vec,
                                 field::GF2_128 delta, field::GF2_128* B) {
    std::vector<field::GF2_128> inverse_input_vole_value(NUM_SBOXES);
    std::vector<field::GF2_128> inverse_output_vole_value(NUM_SBOXES);
    std::vector<field::GF2_128> muti_gate_vole_value(3);

    hash_forward_128_verifier(q, q_vec, delta, inverse_input_vole_value.data(),
                              muti_gate_vole_value.data());
    hash_backword_128_verifier(q, inverse_output_vole_value.data(),
                               muti_gate_vole_value.data() + 2);
    B[0] = muti_gate_vole_value[0] * muti_gate_vole_value[1] -
           muti_gate_vole_value[2] * delta;
    for (int i = 0; i < NUM_SBOXES; i++) {
        B[i + 1] = inverse_input_vole_value[i] * inverse_output_vole_value[i] -
                   delta * delta;
    }
}

void convert_vec_to_field(uint8_t** vec, field::GF2_128* field_vec,
                          const unsigned int ell, const unsigned int lambda) {
    for (unsigned int row = 0; row < ell; row++) {
        uint8_t new_row[lambda / 8] = {0};
        for (unsigned int column = 0; column < 128; ++column) {
            ptr_set_bit(new_row, ptr_get_bit(vec[column], row), column);
        }
        field_vec[row / 8 * 8 + (row % 8)].from_bytes(new_row);
    }
}

void gen_combined_field_vec(field::GF2_128* field_vec,
                            field::GF2_128* combined_field_vec,
                            const unsigned int ell) {
    for (unsigned int i = 0; i < ell / 128; i += 1) {
        combined_field_vec[i] = combine_bf128_vec(field_vec);
        field_vec += 128ULL;
    }
}

// Rain3:
// void path_prove(const uint8_t* witness, field::GF2_128* v,
//                 field::GF2_128* v_vec, const std::vector<uint8_t>& in,
//                 field::GF2_128* A_0, field::GF2_128* A_1,
//                 unsigned int hash_times) {

//     for (int i = 0; i < hash_times; i++) {
//         hash_constrain_128_prover(v + 6 * i, v_vec + (6 * i) * 128,
//                                   witness + (12 * i) * sizeof(uint64_t),
//                                   A_0 + (4 * i), A_1 + (4 * i));
//     }
// }

// void path_verify(field::GF2_128* q, field::GF2_128* q_vec, field::GF2_128
// delta,
//                  const std::vector<uint8_t>& in, field::GF2_128* B,
//                  unsigned int hash_times) {

//     for (int i = 0; i < hash_times; i++) {
//         hash_constrain_128_verifier(q + 6 * i, q_vec + (6 * i) * 128,
//                                     delta, B + (4 * i));
//     }
// }

void path_prove(const uint8_t* witness, field::GF2_128* v,
                field::GF2_128* v_vec, const std::vector<uint8_t>& in,
                field::GF2_128* A_0, field::GF2_128* A_1,
                unsigned int hash_times) {
    for (int i = 0; i < hash_times; i++) {
        hash_constrain_128_prover(v + 7 * i, v_vec + (7 * i) * 128,
                                  witness + (14 * i) * sizeof(uint64_t),
                                  A_0 + (5 * i), A_1 + (5 * i));
    }
}

void path_verify(field::GF2_128* q, field::GF2_128* q_vec, field::GF2_128 delta,
                 const std::vector<uint8_t>& in, field::GF2_128* B,
                 unsigned int hash_times) {
    for (int i = 0; i < hash_times; i++) {
        hash_constrain_128_verifier(q + 7 * i, q_vec + (7 * i) * 128, delta,
                                    B + (5 * i));
    }
}