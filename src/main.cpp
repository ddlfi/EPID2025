#include <chrono>
#include <cstdint>
#include <iostream>

#include "signature.h"

template <typename GF2>
void test_case(const std::vector<uint8_t>& msg, int lambda, int k0, int k1,
               int tau0, int tau1, int tau, int m, int k, int hamming_weight) {
    std::cout << "test params:" << std::endl;
    std::cout << "lambda: " << lambda << " k0: " << k0 << " k1: " << k1
              << " tau0: " << tau0 << " tau1: " << tau1 << " tau: " << tau
              << " m: " << m << " k: " << k
              << " hamming_weight: " << hamming_weight << std::endl;

    Signature<GF2> signer(lambda, k0, k1, tau0, tau1, tau, m, k,
                          hamming_weight);

    signature_t<GF2> sig;
    signer.sign(msg, &sig);

    if (signer.verify(msg, &sig)) {
        std::cout << "pass" << std::endl;
    } else {
        std::cout << "fail" << std::endl;
    }
    std::cout<<std::endl;
}
int main() {
    gen_field_base_128(field_base_128);
    gen_field_base_256(field_base_256);
    std::vector<uint8_t> msg = {0x11, 0x22};

    test_case<field::GF2_256>(msg, 256, 8, 8, 0, 32, 32, 2432, 1216, 258);
    test_case<field::GF2_256>(msg, 256, 12, 11, 14, 8, 22, 2432, 1216, 258);

    test_case<field::GF2_128>(msg, 128, 8, 8, 0, 16, 16, 1268, 634, 133);
    test_case<field::GF2_128>(msg, 128, 12, 11, 7, 4, 11, 1268, 634, 133);
}