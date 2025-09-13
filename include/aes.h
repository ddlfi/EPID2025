/*
 *  SPDX-License-Identifier: MIT
 */

#ifndef FAEST_AES_H
#define FAEST_AES_H

#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include "instances.h"
#include "macros.h"

FAEST_BEGIN_C_DECL

#define AES_MAX_ROUNDS 14

typedef uint8_t aes_word_t[4];
// round key with 4 (AES) up to 8 (Rijndael-256) units
// TODO: aes_round_key_t[8] should this be 8 ?
typedef aes_word_t aes_round_key_t[8];

// # of rows
#define AES_NR 4

// block with 4 (AES) up to 8 (Rijndael-256) units
typedef aes_word_t aes_block_t[8];

typedef struct {
    aes_round_key_t round_keys[AES_MAX_ROUNDS + 1];
} aes_round_keys_t;

int aes128_init_round_keys(aes_round_keys_t* round_key, const uint8_t* key);

int aes128_encrypt_block(const aes_round_keys_t* key, const uint8_t* plaintext,
                         uint8_t* ciphertext);

void aes_increment_iv(uint8_t* iv);

void aes_extend_witness(const uint8_t* key, const uint8_t* in, uint8_t* w,
                        bool origin_key_flag, bool round_key_flag,
                        bool out_flag);

int expand_key(aes_round_keys_t* round_keys, const uint8_t* key,
                   unsigned int key_words, unsigned int block_words,
                   unsigned int num_rounds);

void prg(const uint8_t* key, const uint8_t* iv, uint8_t* out, unsigned int bits,
         size_t outlen);
FAEST_END_C_DECL

#endif