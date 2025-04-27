#include "field.h"
#include "rain_4.h"
#include "utils.h"

void commit_forward_128(const uint8_t* witness, bool flag,
                        field::GF2_128* bf_y);

void commit_forward_128_prover(field::GF2_128* v, field::GF2_128* v_vec,
                               bool flag, field::GF2_128* bf_y);

void commit_forward_128_verifier(field::GF2_128* q, field::GF2_128* q_vec,
                                 field::GF2_128 delta, bool flag,
                                 field::GF2_128* bf_y);

void commit_backword_128(const uint8_t* witness, bool flag,
                         field::GF2_128* bf_y);

void commit_backword_128_prover(field::GF2_128* v, bool flag,
                                field::GF2_128* bf_y);

void commit_backword_128_verifier(field::GF2_128* q, bool flag,
                                  field::GF2_128* bf_y);

void commit_constrain_128_prover(field::GF2_128* v, field::GF2_128* v_vec,
                                 const uint8_t* witness, bool flag,
                                 field::GF2_128* A_0, field::GF2_128* A_1);

void commit_constrain_128_verifier(field::GF2_128* q, field::GF2_128* q_vec,
                                   field::GF2_128 delta, bool flag,
                                   field::GF2_128* B);

void baps_commit_prove(const uint8_t* witness, field::GF2_128* v,
                       field::GF2_128* v_vec, field::GF2_128* A_0,
                       field::GF2_128* A_1);
                       
void baps_commit_verify(field::GF2_128* q, field::GF2_128* q_vec,
                        field::GF2_128 delta, field::GF2_128* B);