#include <chrono>
#include <cstdint>
#include <iostream>

#include "signature.h"

void test_case(unsigned int member_num, unsigned int signer_index,
               const std::vector<uint8_t>& msg) {
    Signature signer(member_num);
    std::cout<<"Start test case : "<<member_num<<" members :"<<std::endl;
    std::cout<<std::endl;
    // std::cout << "         sign stage :" << std::endl;
    auto time_1 = std::chrono::high_resolution_clock::now();
    signer.gen_tree();
    auto time_2 = std::chrono::high_resolution_clock::now();
    auto tree_time = time_2 - time_1;
    // std::cout << "merkle tree generation time is : "
    //           << std::chrono::duration<double, std::milli>(tree_time).count()
    //           << " ms" << std::endl;

    signature_t sig;
    signer.sign(signer_index, msg, &sig);
    auto time_3 = std::chrono::high_resolution_clock::now();

    auto sign_time = time_3 - time_1;
    std::cout << "total sign time is : "
              << std::chrono::duration<double, std::milli>(sign_time).count()
              << " ms" << std::endl;

    // std::cout << "         verify stage :" << std::endl;
    // std::cout << "merkle tree generation time is : "
    //           << std::chrono::duration<double, std::milli>(tree_time).count()
    //           << " ms" << std::endl;
    std::chrono::time_point<std::chrono::high_resolution_clock> time_4;
    if (signer.verify(msg, &sig)) {
        time_4 = std::chrono::high_resolution_clock::now();
        std::cout << "total verify time is : "
                  << std::chrono::duration<double, std::milli>(time_4 - time_3 + time_2 - time_1)
                         .count()
                  << " ms" << std::endl;
        std::cout << "test case : " << member_num << " members pass!"
                  << std::endl;
    } else {
        time_4 = std::chrono::high_resolution_clock::now();
        std::cout << "total sign time is : "
                  << std::chrono::duration<double, std::milli>(time_4 - time_3)
                         .count()
                  << " ms" << std::endl;
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