#ifndef FAEST_UTILS_H
#define FAEST_UTILS_H

#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include "compat.h"
#include "macros.h"

FAEST_BEGIN_C_DECL

static inline void xor_u8_array(const uint8_t* a, const uint8_t* b, uint8_t* out, size_t len) {
  for (size_t i = 0; i < len; i++) {
    out[i] = a[i] ^ b[i];
  }
}

static inline void masked_xor_u8_array(const uint8_t* a, const uint8_t* b, uint8_t* out,
                                       uint8_t mask_bit, size_t len) {
  uint8_t mask = -(mask_bit & 1);
  for (size_t i = 0; i < len; i++) {
    out[i] = a[i] ^ (b[i] & mask);
  }
}

static inline void setbit(uint8_t* value, uint8_t position, uint8_t value_to_set) {
    assert(position < 8);
    const uint8_t mask = 1 << position;
    *value = (*value & ~mask) | ((value_to_set ? 1 : 0) << position);
}

static inline uint8_t parity(uint8_t x) {
    x ^= x >> 4;    
    x ^= x >> 2;    
    x ^= x >> 1;    
    return x & 1;
}

#define get_bit(value, index) (((value) >> (index)) & 1)
#define set_bit(value, index) ((value) << (index))
#define ptr_get_bit(value, index) (((value)[(index) / 8] >> ((index) % 8)) & 1)
#define ptr_set_bit(dst, value, index)                                                             \
  do {                                                                                             \
    (dst)[(index) / 8] |= (value) << ((index) % 8);                                                \
  } while (0)

FAEST_END_C_DECL

#endif
