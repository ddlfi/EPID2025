/*
 * Parallel optimized VOLE operations header
 * SPDX-License-Identifier: MIT
 */

#ifndef FAEST_VOLE_PARALLEL_H
#define FAEST_VOLE_PARALLEL_H

#include "vole.h"

#ifdef __cplusplus
extern "C" {
#endif

// Parallel optimized versions of VOLE operations
void vole_commit_parallel(const uint8_t* rootKey, const uint8_t* iv, unsigned int ellhat,
                         const paramset_t* params, uint8_t* hcom, vec_com_t* vecCom,
                         uint8_t* c, uint8_t* u, uint8_t** v);

void vole_reconstruct_parallel(const uint8_t* iv, const uint8_t* chall,
                              const uint8_t* const* pdec, const uint8_t* const* com_j,
                              uint8_t* hcom, uint8_t** q, unsigned int ellhat,
                              const paramset_t* params);

#ifdef __cplusplus
}
#endif

#endif // FAEST_VOLE_PARALLEL_H