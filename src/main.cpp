#include <cstdint>
#include <iostream>

#include "signature.h"

void test_case(unsigned int member_num, unsigned int signer_index,
               const std::vector<uint8_t>& msg) {
    Signature signer(member_num);
    std::cout << "Start test case : " << member_num << " members:" << std::endl;
    std::cout << std::endl;
    signature_t sig;
    signer.sign(signer_index, msg, &sig);

    if (signer.verify(msg, &sig)) {
        std::cout << "test case : " << member_num << " members pass!"
                  << std::endl;
    } else {
        std::cout << "test case : " << member_num << " members fail!"
                  << std::endl;
    }
    std::cout << std::endl;
}

int main() {
    gen_field_base_128(field_base_128);
    std::vector<uint8_t> msg = {0x11, 0x22};

    test_case(8U, 0U, msg);
    test_case(16U, 8U, msg);
    test_case(32U, 16U, msg);
    test_case(64U, 32U, msg);
    test_case(128U, 64U, msg);
    test_case(256U, 128U, msg);
    test_case(512U, 256U, msg);
    test_case(1024U, 512U, msg);
    test_case(2048U, 1024U, msg);
    test_case(4096U, 2048U, msg);
}