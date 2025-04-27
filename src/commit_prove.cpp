#include "commit_prove.h"

void commit_forward_128(const uint8_t* witness, bool flag,
                        field::GF2_128* bf_y) {
    field::GF2_128 r, m_1;
    r.from_bytes(witness);
    m_1.from_bytes(witness + 16UL);
    witness += 16UL * 2;

    field::GF2_128 inverse_input_value;
    inverse_input_value = r;
    inverse_input_value += m_1;
    inverse_input_value += roundconst[0];
    bf_y[0] = inverse_input_value;

    for (int i = 1; i < NUM_SBOXES; i++) {
        field::GF2_128 tmp;
        tmp.from_bytes(witness);
        witness += 16UL;
        tmp = tmp.multiply_with_transposed_GF2_matrix(matrix_transposed[i - 1]);
        tmp += m_1;
        tmp += roundconst[i];
        bf_y[i] = tmp;
    }
    if (flag) {
        field::GF2_128 out1, m_2;
        out1.from_bytes(witness);
        m_2.from_bytes(witness + 16UL);
        witness += 16UL * 2;

        inverse_input_value = out1;
        inverse_input_value += m_2;
        inverse_input_value += roundconst[0];
        bf_y[4] = inverse_input_value;

        for (int i = 1; i < NUM_SBOXES; i++) {
            field::GF2_128 tmp;
            tmp.from_bytes(witness);
            witness += 16UL;
            tmp = tmp.multiply_with_transposed_GF2_matrix(
                matrix_transposed[i - 1]);
            tmp += m_2;
            tmp += roundconst[i];
            bf_y[4 + i] = tmp;
        }
    }
}

void commit_forward_128_prover(field::GF2_128* v, field::GF2_128* v_vec,
                               bool flag, field::GF2_128* bf_y) {
    bf_y[0] = v[0] + v[1];
    for (int i = 1; i < NUM_SBOXES; i++) {
        bf_y[i] = gf128_vec_muti_with_transposed_GF2_matrix(
            v_vec + 128 * (i + 1), matrix_transposed[i - 1]);
        bf_y[i] += v[1];
    }
    if (flag) {
        bf_y[4] = v[5] + v[6];
        for (int i = 1; i < NUM_SBOXES; i++) {
            bf_y[i + 4] = gf128_vec_muti_with_transposed_GF2_matrix(
                v_vec + 128 * (i + 6), matrix_transposed[i - 1]);
            bf_y[i + 4] += v[6];
        }
    }
}

void commit_forward_128_verifier(field::GF2_128* q, field::GF2_128* q_vec,
                                 field::GF2_128 delta, bool flag,
                                 field::GF2_128* bf_y) {
    bf_y[0] = q[0] + q[1] + delta * roundconst[0];
    for (int i = 1; i < NUM_SBOXES; i++) {
        bf_y[i] = gf128_vec_muti_with_transposed_GF2_matrix(
            q_vec + 128 * (i + 1), matrix_transposed[i - 1]);
        bf_y[i] += delta * roundconst[i] + q[1];
    }
    if (flag) {
        bf_y[4] = q[5] + q[6] + delta * roundconst[0];
        for (int i = 1; i < NUM_SBOXES; i++) {
            bf_y[i + 4] = gf128_vec_muti_with_transposed_GF2_matrix(
                q_vec + 128 * (i + 6), matrix_transposed[i - 1]);
            bf_y[i + 4] += delta * roundconst[i] + q[6];
        }
    }
}

void commit_backword_128(const uint8_t* witness, bool flag,
                         field::GF2_128* bf_y) {
    field::GF2_128 r, m_1, out1;
    r.from_bytes(witness);
    m_1.from_bytes(witness + 16UL);
    witness += 16UL * 2;

    field::GF2_128 tmp;
    for (int i = 0; i < NUM_SBOXES - 1; i++) {
        tmp.from_bytes(witness);
        witness += 16UL;
        bf_y[i] = tmp;
    }

    out1.from_bytes(witness);
    witness += 16UL;
    tmp = out1 + r + m_1;
    bf_y[3] = tmp;

    if (flag) {
        field::GF2_128 m_2;
        m_2.from_bytes(witness);
        witness += 16UL;

        for (int i = 0; i < NUM_SBOXES - 1; i++) {
            tmp.from_bytes(witness);
            witness += 16UL;
            bf_y[i + 4] = tmp;
        }
        tmp.from_bytes(witness);
        witness += 16UL;
        tmp += out1 + m_2;
        bf_y[7] = tmp;
    }
}

void commit_backword_128_prover(field::GF2_128* v, bool flag,
                                field::GF2_128* bf_y) {
    for (int i = 0; i < NUM_SBOXES - 1; i++) {
        bf_y[i] = v[i + 2];
    }
    bf_y[3] = v[0] + v[1] + v[5];
    if (flag) {
        for (int i = 0; i < NUM_SBOXES - 1; i++) {
            bf_y[i + 4] = v[i + 7];
        }
        bf_y[7] = v[5] + v[6] + v[10];
    }
}

void commit_backword_128_verifier(field::GF2_128* q, bool flag,
                                  field::GF2_128* bf_y) {
    for (int i = 0; i < NUM_SBOXES - 1; i++) {
        bf_y[i] = q[i + 2];
    }
    bf_y[3] = q[0] + q[1] + q[5];
    if (flag) {
        for (int i = 0; i < NUM_SBOXES - 1; i++) {
            bf_y[i + 4] = q[i + 7];
        }
        bf_y[7] = q[5] + q[6] + q[10];
    }
}

void commit_constrain_128_prover(field::GF2_128* v, field::GF2_128* v_vec,
                                 const uint8_t* witness, bool flag,
                                 field::GF2_128* A_0, field::GF2_128* A_1) {
    uint8_t NUM_INVERSE = NUM_SBOXES + flag * NUM_SBOXES;
    std::vector<field::GF2_128> inverse_input_real_value(NUM_INVERSE);
    std::vector<field::GF2_128> inverse_output_real_value(NUM_INVERSE);
    std::vector<field::GF2_128> inverse_input_vole_value(NUM_INVERSE);
    std::vector<field::GF2_128> inverse_output_vole_value(NUM_INVERSE);

    commit_forward_128(witness, flag, inverse_input_real_value.data());
    commit_forward_128_prover(v, v_vec, flag, inverse_input_vole_value.data());
    commit_backword_128(witness, flag, inverse_output_real_value.data());
    commit_backword_128_prover(v, flag, inverse_output_vole_value.data());

    for (int i = 0; i < NUM_INVERSE; i++) {
        A_0[i] = inverse_input_vole_value[i] * inverse_output_vole_value[i];
        A_1[i] = inverse_input_vole_value[i] * inverse_output_real_value[i] +
                 inverse_output_vole_value[i] * inverse_input_real_value[i];
    }
}

void commit_constrain_128_verifier(field::GF2_128* q, field::GF2_128* q_vec,
                                   field::GF2_128 delta, bool flag,
                                   field::GF2_128* B) {
    uint8_t NUM_INVERSE = NUM_SBOXES + flag * NUM_SBOXES;
    std::vector<field::GF2_128> inverse_input_vole_value(NUM_INVERSE);
    std::vector<field::GF2_128> inverse_output_vole_value(NUM_INVERSE);

    commit_forward_128_verifier(q, q_vec, delta, flag,
                                inverse_input_vole_value.data());
    commit_backword_128_verifier(q, flag, inverse_output_vole_value.data());
    for (int i = 0; i < NUM_INVERSE; i++) {
        B[i] = inverse_input_vole_value[i] * inverse_output_vole_value[i] -
               delta * delta;
    }
}

void baps_commit_prove(const uint8_t* witness, field::GF2_128* v,
                       field::GF2_128* v_vec, field::GF2_128* A_0,
                       field::GF2_128* A_1) {
    commit_constrain_128_prover(v, v_vec, witness, 0, A_0, A_1);
    commit_constrain_128_prover(v + 6, v_vec + 6 * 128, witness + 6 * 16UL, 1,
                                A_0 + 4, A_1 + 4);
}

void baps_commit_verify(field::GF2_128* q, field::GF2_128* q_vec,
                        field::GF2_128 delta, field::GF2_128* B) {
    commit_constrain_128_verifier(q, q_vec, delta, 0, B);
    commit_constrain_128_verifier(q + 6, q_vec + 6 * 128, delta, 1, B + 4);
}