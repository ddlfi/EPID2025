#include <chrono>
#include <cstdint>
#include <iostream>

#include "epid.h"

void test_case(const std::vector<uint8_t>& msg, int lambda, int k0, int k1,
               int tau0, int tau1, int tau, int member_num,int test_srl_size) {
    std::cout << "test params:" << std::endl;
    std::cout << "lambda: " << lambda << " k0: " << k0 << " k1: " << k1
              << " tau0: " << tau0 << " tau1: " << tau1 << " tau: " << tau
              << std::endl;
    
    std::cout << "Generating simulation environment with " << member_num << " members and SRL size " << test_srl_size << "..." << std::endl;
    // Initialize simulation environment, generate group members and SRL of specified size
    epid signer(lambda, k0, k1, tau0, tau1, tau, member_num,test_srl_size);
    signature_t sig;
    
    std::cout << "Starting signature generation..." << std::endl;
    signer.epid_sign(msg, &sig, 1);

    std::cout << "Starting signature verification..." << std::endl;
    if (signer.epid_verify(msg, &sig)) {
        std::cout << "pass" << std::endl;
    } else {
        std::cout << "fail" << std::endl;
    }
    std::cout << std::endl;
}
int main() {
    std::vector<uint8_t> msg = {0x11, 0x22};

    test_case(msg, 128, 8, 8, 0, 16, 16, 1024*1024, 1000);
}