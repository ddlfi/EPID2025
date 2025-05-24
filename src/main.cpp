#include <chrono>
#include <cstdint>
#include <iostream>

#include "signature.h"

void test_case(const std::vector<uint8_t>& msg) {
    Signature signer;

    signature_t sig;
    signer.sign(msg, &sig);

    if (signer.verify(msg, &sig)) {
        std::cout << "pass" << std::endl;
    } else {
        std::cout << "fail" << std::endl;
    }
}
int main() {
    gen_field_base_128(field_base_128);
    std::vector<uint8_t> msg = {0x11, 0x22};
    test_case(msg);
}