/*
 * Parallel optimized version of VOLE operations
 * SPDX-License-Identifier: MIT
 */

#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#include "aes.h"
#include "random_oracle.h"
#include "utils.h"
#include "vole.h"
#include "vc.h"
#include "vole_parallel.h"

#ifndef MIN
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#endif

// Memory-optimized version of ConvertToVole that uses less memory
static void ConvertToVole_memory_optimized(const uint8_t* iv, const uint8_t* sd, bool sd0_bot,
                                          unsigned int lambda, unsigned int depth,
                                          unsigned int outLenBytes, uint8_t* u, uint8_t* v) {
    const unsigned int num_instances = 1 << depth;
    const unsigned int lambda_bytes = lambda / 8;

    // Allocate smaller buffer - only what we need at once
    // Instead of 2 * num_instances, use smaller chunks
    const unsigned int chunk_size = 256; // Process in chunks of 256 instances
    uint8_t* r = calloc(2 * chunk_size, outLenBytes);
    uint8_t* temp_r = calloc(num_instances, outLenBytes); // Temporary storage for one row

#define R_CHUNK(row, column) \
    (r + (((row) % 2) * chunk_size + (column)) * outLenBytes)
#define TEMP_R(column) \
    (temp_r + (column) * outLenBytes)
#define V(idx) (v + (idx) * outLenBytes)

    // Step: 2
    if (!sd0_bot) {
        prg(sd, iv, TEMP_R(0), lambda, outLenBytes);
    }
    
    // Step: 3..4 - Process in chunks to save memory
    for (unsigned int chunk = 0; chunk < (num_instances + chunk_size - 1) / chunk_size; chunk++) {
        unsigned int start_idx = chunk * chunk_size;
        unsigned int end_idx = MIN(start_idx + chunk_size, num_instances);
        
        for (unsigned int i = start_idx; i < end_idx; i++) {
            if (i == 0 && !sd0_bot) {
                // Already computed above
                continue;
            }
            prg(sd + (lambda_bytes * i), iv, TEMP_R(i), lambda, outLenBytes);
        }
    }
    
    // Step: 5..9 - Process level by level with memory-efficient approach
    memset(v, 0, depth * outLenBytes);
    
    // Copy initial row from temp_r to processing buffer
    for (unsigned int i = 0; i < num_instances; i++) {
        memcpy(R_CHUNK(0, i % chunk_size), TEMP_R(i), outLenBytes);
        
        // Process when chunk is full or at end
        if ((i + 1) % chunk_size == 0 || i == num_instances - 1) {
            unsigned int chunk_start = (i / chunk_size) * chunk_size;
            unsigned int chunk_end = MIN(chunk_start + chunk_size, num_instances);
            
            // Process this chunk for all depth levels
            for (unsigned int j = 0; j < depth; j++) {
                unsigned int level_instances = num_instances >> (j + 1);
                unsigned int chunk_level_start = chunk_start >> (j + 1);
                unsigned int chunk_level_end = MIN(chunk_end >> (j + 1), level_instances);
                
                for (unsigned int k = chunk_level_start; k < chunk_level_end; k++) {
                    unsigned int global_k = k - chunk_level_start;
                    if (2 * k + 1 < chunk_end) {
                        xor_u8_array(V(j), R_CHUNK(j, 2 * global_k + 1), V(j), outLenBytes);
                        xor_u8_array(R_CHUNK(j, 2 * global_k), R_CHUNK(j, 2 * global_k + 1), 
                                   R_CHUNK(j + 1, global_k), outLenBytes);
                    }
                }
            }
        }
    }
    
    // Step: 10
    if (!sd0_bot && u != NULL) {
        memcpy(u, R_CHUNK(depth, 0), outLenBytes);
    }
    
    free(r);
    free(temp_r);

#undef R_CHUNK
#undef TEMP_R
#undef V
}

void vole_commit_parallel(const uint8_t* rootKey, const uint8_t* iv, unsigned int ellhat,
                         const paramset_t* params, uint8_t* hcom, vec_com_t* vecCom,
                         uint8_t* c, uint8_t* u, uint8_t** v) {
    unsigned int lambda = params->lambda;
    unsigned int lambda_bytes = lambda / 8;
    unsigned int ellhat_bytes = ellhat / 8;
    unsigned int tau = params->tau;
    unsigned int tau0 = params->tau0;
    unsigned int k0 = params->k0;
    unsigned int k1 = params->k1;

    // Calculate memory requirements for ConvertToVole
    unsigned int max_depth = MAX(k0, k1);
    unsigned int max_num_instances = 1 << max_depth;
    size_t memory_per_convertovole = 2ULL * max_num_instances * ellhat_bytes;
    size_t total_parallel_memory = tau * memory_per_convertovole;
    
    // Check available memory
    FILE* meminfo = fopen("/proc/meminfo", "r");
    size_t available_memory = 0;
    if (meminfo) {
        char line[256];
        while (fgets(line, sizeof(line), meminfo)) {
            if (strncmp(line, "MemAvailable:", 13) == 0) {
                unsigned long available_kb;
                sscanf(line, "MemAvailable: %lu kB", &available_kb);
                available_memory = available_kb * 1024ULL; // Convert to bytes
                break;
            }
        }
        fclose(meminfo);
    }
    
    // Decision: use parallel if we have enough memory (with 20% safety margin)
    bool use_parallel = (available_memory > 0) && (total_parallel_memory < available_memory * 0.8);
    
    printf("    [VOLE] Memory analysis:\n");
    printf("    [VOLE]   Required per ConvertToVole: %.2f GB\n", 
           memory_per_convertovole / (1024.0 * 1024.0 * 1024.0));
    printf("    [VOLE]   Total parallel requirement: %.2f GB\n", 
           total_parallel_memory / (1024.0 * 1024.0 * 1024.0));
    printf("    [VOLE]   Available memory: %.2f GB\n", 
           available_memory / (1024.0 * 1024.0 * 1024.0));
    printf("    [VOLE]   Strategy: %s\n", use_parallel ? "PARALLEL" : "SERIAL");
    fflush(stdout);

    uint8_t* ui = malloc(tau * ellhat_bytes);

    // Step 1: Key expansion (serial)
    uint8_t* expanded_keys = malloc(tau * lambda_bytes);
    prg(rootKey, iv, expanded_keys, lambda, lambda_bytes * tau);

    // Precompute v indices
    unsigned int* v_indices = malloc(tau * sizeof(unsigned int));
    unsigned int v_idx = 0;
    for (unsigned int i = 0; i < tau; i++) {
        v_indices[i] = v_idx;
        unsigned int depth = i < tau0 ? k0 : k1;
        v_idx += depth;
    }

    if (use_parallel) {
        // Parallel version: both vector_commitment and ConvertToVole in parallel
        printf("    [VOLE] Using full parallel implementation\n");
        fflush(stdout);
        
        #pragma omp parallel for schedule(dynamic)
        for (unsigned int i = 0; i < tau; i++) {
            unsigned int depth = i < tau0 ? k0 : k1;
            
            // Step 5: vector_commitment
            vector_commitment(expanded_keys + i * lambda_bytes, iv, params, lambda,
                              &vecCom[i], depth);
            
            // Step 6: ConvertToVole (original version, high memory but faster)
            ConvertToVole(iv, vecCom[i].sd, false, lambda, depth, ellhat_bytes,
                          ui + i * ellhat_bytes, v[v_indices[i]]);
        }
    } else {
        // Memory-optimized version: serial outer loop with optimized ConvertToVole
        printf("    [VOLE] Using memory-optimized implementation\n");
        fflush(stdout);
        
        // Phase 1: Parallel vector_commitment (low memory usage)
        #pragma omp parallel for schedule(dynamic)
        for (unsigned int i = 0; i < tau; i++) {
            unsigned int depth = i < tau0 ? k0 : k1;
            vector_commitment(expanded_keys + i * lambda_bytes, iv, params, lambda,
                              &vecCom[i], depth);
        }
        
        // Phase 2: Serial ConvertToVole_memory_optimized (lower memory usage)
        for (unsigned int i = 0; i < tau; i++) {
            unsigned int depth = i < tau0 ? k0 : k1;
            ConvertToVole_memory_optimized(iv, vecCom[i].sd, false, lambda, depth, ellhat_bytes,
                                         ui + i * ellhat_bytes, v[v_indices[i]]);
        }
    }

    free(expanded_keys);
    free(v_indices);

    // Step 9-11: XOR operations (serial)
    memcpy(u, ui, ellhat_bytes);
    for (unsigned int i = 1; i < tau; i++) {
        xor_u8_array(u, ui + i * ellhat_bytes, c + (i - 1) * ellhat_bytes,
                     ellhat_bytes);
    }
    free(ui);

    // Step 12: H1 hashing (serial, order-dependent)
    H1_context_t h1_ctx;
    H1_init(&h1_ctx, lambda);
    for (unsigned int i = 0; i < tau; i++) {
        H1_update(&h1_ctx, vecCom[i].h, lambda_bytes * 2);
    }
    H1_final(&h1_ctx, hcom, lambda_bytes * 2);
}

void vole_reconstruct_parallel(const uint8_t* iv, const uint8_t* chall,
                              const uint8_t* const* pdec, const uint8_t* const* com_j,
                              uint8_t* hcom, uint8_t** q, unsigned int ellhat,
                              const paramset_t* params) {
    unsigned int lambda = params->lambda;
    unsigned int lambda_bytes = lambda / 8;
    unsigned int ellhat_bytes = (ellhat + 7) / 8;
    unsigned int tau = params->tau;
    unsigned int tau0 = params->tau0;
    unsigned int tau1 = params->tau1;
    unsigned int k0 = params->k0;
    unsigned int k1 = params->k1;

    // Same memory analysis as vole_commit_parallel
    unsigned int max_depth = MAX(k0, k1);
    unsigned int max_num_instances = 1 << max_depth;
    size_t memory_per_convertovole = 2ULL * max_num_instances * ellhat_bytes;
    size_t total_parallel_memory = tau * memory_per_convertovole;
    
    // Check available memory
    FILE* meminfo = fopen("/proc/meminfo", "r");
    size_t available_memory = 0;
    if (meminfo) {
        char line[256];
        while (fgets(line, sizeof(line), meminfo)) {
            if (strncmp(line, "MemAvailable:", 13) == 0) {
                unsigned long available_kb;
                sscanf(line, "MemAvailable: %lu kB", &available_kb);
                available_memory = available_kb * 1024ULL;
                break;
            }
        }
        fclose(meminfo);
    }
    
    bool use_parallel = (available_memory > 0) && (total_parallel_memory < available_memory * 0.8);
    
    printf("    [VOLE_RECONSTRUCT] Strategy: %s\n", use_parallel ? "PARALLEL" : "SERIAL");
    fflush(stdout);

    if (use_parallel) {
        // Full parallel version
        uint8_t** sd_per_thread = malloc(tau * sizeof(uint8_t*));
        vec_com_rec_t* vecComRec_per_thread = malloc(tau * sizeof(vec_com_rec_t));
        
        for (unsigned int i = 0; i < tau; i++) {
            sd_per_thread[i] = malloc((1 << MAX(k0, k1)) * lambda_bytes);
            memset(sd_per_thread[i], 0, lambda_bytes);
            
            vecComRec_per_thread[i].h = malloc(lambda_bytes * 2);
            vecComRec_per_thread[i].k = calloc(getBinaryTreeNodeCount(MAX(k0, k1)), lambda_bytes);
            vecComRec_per_thread[i].com = malloc((1 << MAX(k0, k1)) * lambda_bytes * 2);
            vecComRec_per_thread[i].s = malloc((1 << MAX(k0, k1)) * lambda_bytes);
        }

        unsigned int* q_indices = malloc(tau * sizeof(unsigned int));
        unsigned int q_idx = 0;
        for (unsigned int i = 0; i < tau; i++) {
            q_indices[i] = q_idx;
            unsigned int depth = i < tau0 ? k0 : k1;
            q_idx += depth;
        }

        // Parallel processing
        #pragma omp parallel for schedule(dynamic)
        for (unsigned int i = 0; i < tau; i++) {
            unsigned int depth = i < tau0 ? k0 : k1;
            unsigned int N = 1 << depth;

            uint8_t chalout[MAX_DEPTH];
            ChalDec(chall, i, k0, tau0, k1, tau1, chalout);
            unsigned int idx = NumRec(depth, chalout);

            vector_reconstruction(iv, pdec[i], com_j[i], chalout, lambda, depth,
                                  &vecComRec_per_thread[i]);

            for (unsigned int j = 1; j < N; j++) {
                memcpy(sd_per_thread[i] + j * lambda_bytes,
                       vecComRec_per_thread[i].s + (lambda_bytes * (j ^ idx)), lambda_bytes);
            }

            ConvertToVole(iv, sd_per_thread[i], true, lambda, depth, ellhat_bytes, NULL,
                          q[q_indices[i]]);
        }

        // Serial H1 hashing
        H1_context_t h1_ctx;
        H1_init(&h1_ctx, lambda);
        for (unsigned int i = 0; i < tau; i++) {
            H1_update(&h1_ctx, vecComRec_per_thread[i].h, lambda_bytes * 2);
        }
        H1_final(&h1_ctx, hcom, lambda_bytes * 2);

        // Cleanup
        for (unsigned int i = 0; i < tau; i++) {
            vec_com_rec_clear(&vecComRec_per_thread[i]);
            free(sd_per_thread[i]);
        }
        free(sd_per_thread);
        free(vecComRec_per_thread);
        free(q_indices);
        
    } else {
        // Memory-safe serial version (original algorithm)
        uint8_t* sd = malloc((1 << MAX(k0, k1)) * lambda_bytes);
        memset(sd, 0, lambda_bytes);

        H1_context_t h1_ctx;
        H1_init(&h1_ctx, lambda);

        vec_com_rec_t vecComRec;
        vecComRec.h = malloc(lambda_bytes * 2);
        vecComRec.k = calloc(getBinaryTreeNodeCount(MAX(k0, k1)), lambda_bytes);
        vecComRec.com = malloc((1 << MAX(k0, k1)) * lambda_bytes * 2);
        vecComRec.s = malloc((1 << MAX(k0, k1)) * lambda_bytes);

        unsigned int q_idx = 0;
        for (unsigned int i = 0; i < tau; i++) {
            unsigned int depth = i < tau0 ? k0 : k1;
            unsigned int N = 1 << depth;

            uint8_t chalout[MAX_DEPTH];
            ChalDec(chall, i, k0, tau0, k1, tau1, chalout);
            unsigned int idx = NumRec(depth, chalout);

            vector_reconstruction(iv, pdec[i], com_j[i], chalout, lambda, depth, &vecComRec);

            for (unsigned int j = 1; j < N; j++) {
                memcpy(sd + j * lambda_bytes,
                       vecComRec.s + (lambda_bytes * (j ^ idx)), lambda_bytes);
            }

            ConvertToVole(iv, sd, true, lambda, depth, ellhat_bytes, NULL, q[q_idx]);
            q_idx += depth;

            H1_update(&h1_ctx, vecComRec.h, lambda_bytes * 2);
        }
        
        vec_com_rec_clear(&vecComRec);
        free(sd);
        H1_final(&h1_ctx, hcom, lambda_bytes * 2);
    }
}