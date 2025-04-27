#include "hamming_weight_prove.h"

void gen_witness_255(uint8_t* inputs, uint8_t* value, int n) {
    if (n == 7) {
        std::copy(inputs, inputs + 256 / 8, value);
    } else if (n == 8) {
        std::copy(inputs, inputs + 512 / 8, value);
    }
    int pos = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < (1 << (n - i)) - 1; j++) {
            FA(value + pos / 8, i, j, n);
        }
        HA(value + pos / 8, i, n);
        if (i >= n - 2) {
            pos += 16;
        } else {
            pos += ((1 << (n + 1 - i)) + (1 << (n - i)));
        }
    }
    HA(value + pos / 8, 7, n);
}

void FA(uint8_t* values, int layer, int index, int n) {
    const int pos_1 = std::max(1 << (n + 1 - layer), 8);
    const int pos_2 = std::max(1 << (n - layer), 8);
    uint8_t x1, x2, x3, sumbit, carrybit;
    if (index == 0) {
        x1 = GET_BIT(values[0], 0);
        x2 = GET_BIT(values[0], 1);
        x3 = GET_BIT(values[0], 2);
        sumbit = x1 ^ x2 ^ x3;
        carrybit = ((x1 ^ x2) & (x1 ^ x3)) ^ x1;
    } else {
        x1 = GET_BIT(values[pos_1 / 8 + (index - 1) / 8], (index - 1) % 8);
        x2 = GET_BIT(values[(2 * index + 1) / 8], (2 * index + 1) % 8);
        x3 = GET_BIT(values[(2 * index + 2) / 8], (2 * index + 2) % 8);
        sumbit = x1 ^ x2 ^ x3;
        carrybit = ((x1 ^ x2) & (x1 ^ x3)) ^ x1;
    }
    setbit(values + pos_1 / 8 + index / 8, index % 8, sumbit);
    setbit(values + pos_1 / 8 + pos_2 / 8 + index / 8, index % 8, carrybit);
}

void HA(uint8_t* values, int layer, int n) {
    if (n == 8) return;
    // layer == 7:
    const int pos_1 = std::max(1 << (n + 1 - layer), 8);
    const int pos_2 = std::max(1 << (n - layer), 8);
    uint8_t x1, x2, sumbit, carrybit;
    x1 = GET_BIT(values[pos_1 / 8 - 1], 7);
    if (layer == 7) {
        x2 = GET_BIT(values[pos_1 / 8 - 1], 0);
    } else {
        x2 = GET_BIT(values[pos_1 / 8 + pos_2 / 8 - 1],
                     (1 << (n - layer) - 2) % 8);
    }
    carrybit = x1 & x2;
    sumbit = x1 ^ x2;
    setbit(values + pos_1 / 8 + pos_2 / 8 + pos_2 / 8 - 1, 7, carrybit);
    setbit(values + pos_1 / 8 + pos_2 / 8 - 1, 7, sumbit);
}

void setbit(uint8_t* value, uint8_t position, uint8_t value_to_set) {
    assert(position < 8);
    const uint8_t mask = 1 << position;
    *value = (*value & ~mask) | ((value_to_set ? 1 : 0) << position);
}

void get_mutigate_u(const uint8_t* witness, uint8_t* bf_y, int n) {
    int pos_1 = (1 << (n + 1)), pos_2 = (1 << n);
    int index = 0;
    uint8_t x1, x2, x3;

    for (int i = 0; i < n; i++) {
        x1 = GET_BIT(witness[0], 0);
        x2 = GET_BIT(witness[0], 1);
        x3 = GET_BIT(witness[0], 2);
        setbit(&bf_y[index], 0, x1 ^ x2);
        setbit(&bf_y[index], 1, x1 ^ x3);
        setbit(&bf_y[index], 2, (x1 ^ x2) & (x1 ^ x3));
        index++;
        for (int j = 1; j < (1 << (n - i)) - 1; j++) {
            x1 = GET_BIT(witness[pos_1 / 8 + (j - 1) / 8], (j - 1) % 8);
            x2 = GET_BIT(witness[(2 * j + 1) / 8], (2 * j + 1) % 8);
            x3 = GET_BIT(witness[(2 * j + 2) / 8], (2 * j + 2) % 8);
            setbit(&bf_y[index], 0, x1 ^ x2);
            setbit(&bf_y[index], 1, x1 ^ x3);
            setbit(&bf_y[index], 2, (x1 ^ x2) & (x1 ^ x3));
            index++;
        }
        if (n == 7) {
            x1 = GET_BIT(witness[pos_1 / 8 - 1], 7);
            x2 = GET_BIT(witness[pos_1 / 8 + pos_2 / 8 - 1],
                         ((1 << (n - i)) - 2) % 8);
            setbit(&bf_y[index], 0, x1);
            setbit(&bf_y[index], 1, x2);
            setbit(&bf_y[index], 2, x1 & x2);
            index++;
        }
        witness += pos_1 / 8 + pos_2 / 8;

        if (pos_1 > 8) pos_1 /= 2;
        if (pos_2 > 8) pos_2 /= 2;
    }
    if (n == 7) {
        x1 = GET_BIT(witness[pos_1 / 8 - 1], 7);
        x2 = GET_BIT(witness[pos_1 / 8 - 1], 0);
        setbit(&bf_y[index], 0, x1);
        setbit(&bf_y[index], 1, x2);
        setbit(&bf_y[index], 2, x1 & x2);
    }
    // std::cout << GET_BIT(bf_y[123], 0) << " " << GET_BIT(bf_y[123], 1) << " "
    //           << GET_BIT(bf_y[123], 2) << std::endl;
}

void get_mutigate_v(const field::GF2_128* v, field::GF2_128* bf_y, int n) {
    int pos_1 = (1 << (n + 1)), pos_2 = (1 << n);
    int index = 0;
    field::GF2_128 v1, v2, v3, v_carrybit;
    for (int i = 0; i < n; i++) {
        v1 = v[0];
        v2 = v[1];
        v3 = v[2];
        v_carrybit = v[pos_1 + pos_2];
        bf_y[index++] = v1 + v2;
        bf_y[index++] = v1 + v3;
        bf_y[index++] = v_carrybit + v1;
        for (int j = 1; j < (1 << (n - i)) - 1; j++) {
            v1 = v[pos_1 + j - 1];
            v2 = v[2 * j + 1];
            v3 = v[2 * j + 2];
            v_carrybit = v[pos_1 + pos_2 + j];
            bf_y[index++] = v1 + v2;
            bf_y[index++] = v1 + v3;
            bf_y[index++] = v_carrybit + v1;
        }
        if (n == 7) {
            v1 = v[pos_1 - 1];
            v2 = v[pos_1 + pos_2 - 2];
            v_carrybit = v[pos_1 + pos_2 + pos_2 - 1];
            bf_y[index++] = v1;
            bf_y[index++] = v2;
            bf_y[index++] = v_carrybit;
        }
        v += pos_1 + pos_2;
        if (pos_1 > 8) pos_1 /= 2;
        if (pos_2 > 8) pos_2 /= 2;
    }
    if (n == 7) {
        v1 = v[7];
        v2 = v[0];
        v_carrybit = v[pos_1 + pos_2 + 7];
        bf_y[index++] = v1;
        bf_y[index++] = v2;
        bf_y[index++] = v_carrybit;
    }
}

void get_mutigate_q(const field::GF2_128* q, field::GF2_128* bf_y, int n) {
    int pos_1 = (1 << (n + 1)), pos_2 = (1 << n);
    int index = 0;
    field::GF2_128 q1, q2, q3, q_carrybit;
    for (int i = 0; i < n; i++) {
        q1 = q[0];
        q2 = q[1];
        q3 = q[2];
        q_carrybit = q[pos_1 + pos_2];
        bf_y[index++] = q1 + q2;
        bf_y[index++] = q1 + q3;
        bf_y[index++] = q_carrybit + q1;
        for (int j = 1; j < (1 << (n - i)) - 1; j++) {
            q1 = q[pos_1 + j - 1];
            q2 = q[2 * j + 1];
            q3 = q[2 * j + 2];
            q_carrybit = q[pos_1 + pos_2 + j];
            bf_y[index++] = q1 + q2;
            bf_y[index++] = q1 + q3;
            bf_y[index++] = q_carrybit + q1;
        }
        if (n == 7) {
            q1 = q[pos_1 - 1];
            q2 = q[pos_1 + pos_2 - 2];
            q_carrybit = q[pos_1 + pos_2 + pos_2 - 1];
            bf_y[index++] = q1;
            bf_y[index++] = q2;
            bf_y[index++] = q_carrybit;
        }
        q += pos_1 + pos_2;
        if (pos_1 > 8) pos_1 /= 2;
        if (pos_2 > 8) pos_2 /= 2;
    }
    if (n == 7) {
        q1 = q[7];
        q2 = q[0];
        q_carrybit = q[pos_1 + pos_2 + 7];
        bf_y[index++] = q1;
        bf_y[index++] = q2;
        bf_y[index++] = q_carrybit;
    }
}

void get_constrain_prover(const field::GF2_128* v, const uint8_t* witness,
                          field::GF2_128* A_0, field::GF2_128* A_1,
                          int mutigate_num, int n) {
    std::vector<uint8_t> mutigate_u(mutigate_num);
    std::vector<field::GF2_128> mutigate_v(3 * mutigate_num);
    get_mutigate_u(witness, mutigate_u.data(), n);
    get_mutigate_v(v, mutigate_v.data(), n);
    for (int i = 0; i < mutigate_num; i++) {
        A_0[i] = mutigate_v[3 * i] * mutigate_v[3 * i + 1];
        if (GET_BIT(mutigate_u[i], 0) == 1) A_1[i] += mutigate_v[3 * i + 1];
        if (GET_BIT(mutigate_u[i], 1) == 1) A_1[i] += mutigate_v[3 * i];
        A_1[i] -= mutigate_v[3 * i + 2];
    }
}

void get_constrain_verifier(const field::GF2_128* q, field::GF2_128 delta,
                            field::GF2_128* B, int mutigate_num, int n) {
    std::vector<field::GF2_128> mutigate_q(3 * mutigate_num);
    get_mutigate_q(q, mutigate_q.data(), n);
    for (int i = 0; i < mutigate_num; i++) {
        B[i] = mutigate_q[3 * i] * mutigate_q[3 * i + 1] -
               mutigate_q[3 * i + 2] * delta;
    }
}