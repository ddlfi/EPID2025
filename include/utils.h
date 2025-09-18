#ifndef FAEST_UTILS_H
#define FAEST_UTILS_H

#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include "compat.h"
#include "macros.h"

#ifdef __AVX2__
#include <immintrin.h>
#elif defined(__SSE2__)
#include <emmintrin.h>
#endif

FAEST_BEGIN_C_DECL

static inline void xor_u8_array(const uint8_t* a, const uint8_t* b, uint8_t* out, size_t len) {
#ifdef __AVX2__
    // Process 32 bytes at a time with AVX2
    size_t simd_len = len & ~31;
    for (size_t i = 0; i < simd_len; i += 32) {
        __m256i va = _mm256_loadu_si256((const __m256i*)(a + i));
        __m256i vb = _mm256_loadu_si256((const __m256i*)(b + i));
        __m256i result = _mm256_xor_si256(va, vb);
        _mm256_storeu_si256((__m256i*)(out + i), result);
    }
    // Handle remaining bytes
    for (size_t i = simd_len; i < len; i++) {
        out[i] = a[i] ^ b[i];
    }
#elif defined(__SSE2__)
    // Process 16 bytes at a time with SSE2
    size_t simd_len = len & ~15;
    for (size_t i = 0; i < simd_len; i += 16) {
        __m128i va = _mm_loadu_si128((const __m128i*)(a + i));
        __m128i vb = _mm_loadu_si128((const __m128i*)(b + i));
        __m128i result = _mm_xor_si128(va, vb);
        _mm_storeu_si128((__m128i*)(out + i), result);
    }
    // Handle remaining bytes
    for (size_t i = simd_len; i < len; i++) {
        out[i] = a[i] ^ b[i];
    }
#else
    // Fallback to scalar implementation
    for (size_t i = 0; i < len; i++) {
        out[i] = a[i] ^ b[i];
    }
#endif
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
