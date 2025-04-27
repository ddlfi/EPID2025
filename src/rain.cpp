#include "rain_4.h"

#include <optional>

bool rain(const std::vector<uint8_t> &key_in,
          const std::vector<uint8_t> &plaintext_in,
          std::vector<uint8_t> &ciphertext_out, uint8_t *witness, bool flag) {
    ciphertext_out.resize(BLOCK_SIZE);

    field::GF2_128 key, state;
    key.from_bytes(key_in.data());
    state.from_bytes(plaintext_in.data());

    if (flag && witness) {
        key.to_bytes(witness);
        witness += 16UL;
    }
    // first r-1 rounds
    for (size_t r = 0; r < NUM_SBOXES - 1; r++) {
        state += key;
        state += roundconst[r];
        if (state.is_zero()) return false;
        state = state.inverse();
        // get the w_state
        if (flag && witness) {
            state.to_bytes(witness);
            witness += 16UL;
        }
        // transposed matrix multiplication is faster, so we use that instead
        // standard multiplication can be useful for debugging, set
        // INCLUDE_STANDARD_MATRICES to include them
        //
        // state = state.multiply_with_GF2_matrix(Params::matrix[r]);
        state = state.multiply_with_transposed_GF2_matrix(matrix_transposed[r]);
    }
    // last round
    state += key;
    state += roundconst[NUM_SBOXES - 1];
    if (state.is_zero()) return false;
    state = state.inverse();
    state += key;

    state.to_bytes(ciphertext_out.data());
    if (flag && witness) {
        state.to_bytes(witness);
    }
    return true;
}

bool rain(const field::GF2_128 &key_in, const field::GF2_128 &plaintext_in,
          field::GF2_128 &ciphertext_out, uint8_t *witness, bool flag) {
    field::GF2_128 key, state;
    key = key_in;
    state = plaintext_in;

    // first r-1 rounds
    for (size_t r = 0; r < NUM_SBOXES - 1; r++) {
        state += key;
        state += roundconst[r];
        if (state.is_zero()) return false;
        state = state.inverse();
        // get the w_state
        if (flag && witness) {
            state.to_bytes(witness);
            witness += 16UL;
        }
        // transposed matrix multiplication is faster, so we use that instead
        // standard multiplication can be useful for debugging, set
        // INCLUDE_STANDARD_MATRICES to include them
        //
        // state = state.multiply_with_GF2_matrix(Params::matrix[r]);
        state = state.multiply_with_transposed_GF2_matrix(matrix_transposed[r]);
    }
    // last round
    state += key;
    state += roundconst[NUM_SBOXES - 1];
    if (state.is_zero()) return false;
    state = state.inverse();
    state += key;
    ciphertext_out = state;

    return true;
}
