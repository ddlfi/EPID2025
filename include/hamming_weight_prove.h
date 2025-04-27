#include <iostream>
#include <algorithm>

#include "field.h"
#include "rain_4.h"
#include "utils.h"



void gen_witness_255(uint8_t* inputs, uint8_t* value, int n);

void FA(uint8_t* values, int layer, int index, int n);

void HA(uint8_t* values, int layer, int n);

void setbit(uint8_t* value, uint8_t position, uint8_t value_to_set);

void get_mutigate_u(const uint8_t* witness, uint8_t* bf_y, int n);

void get_mutigate_v(const field::GF2_128* v, field::GF2_128* bf_y, int n);

void get_mutigate_q(const field::GF2_128* q, field::GF2_128* bf_y, int n);

void get_constrain_prover(const field::GF2_128* v, const uint8_t* witness,
                          field::GF2_128* A_0, field::GF2_128* A_1,
                          int mutigate_num,int n);

void get_constrain_verifier(const field::GF2_128* q, field::GF2_128 delta,
                            field::GF2_128* B, int mutigate_num, int n);