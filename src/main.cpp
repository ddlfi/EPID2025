#include <chrono>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>

#include "epid.h"

void test_case(const std::vector<uint8_t>& msg, int lambda, int k0, int k1,
               int tau0, int tau1, int tau, int member_num, int test_srl_size) {
    std::cout << "test params:" << std::endl;
    std::cout << "lambda: " << lambda << " k0: " << k0 << " k1: " << k1
              << " tau0: " << tau0 << " tau1: " << tau1 << " tau: " << tau
              << std::endl;

    std::cout << "Generating simulation environment with " << member_num
              << " members and SRL size " << test_srl_size << "..."
              << std::endl;
    // Initialize simulation environment, generate group members and SRL of
    // specified size
    epid signer(lambda, k0, k1, tau0, tau1, tau, member_num, test_srl_size);
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

void performance_test(const std::vector<uint8_t>& msg, int lambda, int k0,
                      int k1, int tau0, int tau1, int tau, int member_num,
                      int test_srl_size, int iterations = 10) {
    std::cout << "Performance Test - Running " << iterations << " iterations..."
              << std::endl;
    std::cout << "test params:" << std::endl;
    std::cout << "lambda: " << lambda << " k0: " << k0 << " k1: " << k1
              << " tau0: " << tau0 << " tau1: " << tau1 << " tau: " << tau
              << std::endl;

    std::vector<double> sign_times;
    std::vector<double> verify_times;
    epid signer(lambda, k0, k1, tau0, tau1, tau, member_num, test_srl_size);

    for (int i = 0; i < iterations; i++) {
        std::cout << "Iteration " << (i + 1) << "/" << iterations << std::endl;

        signature_t sig;

        auto start_sign = std::chrono::high_resolution_clock::now();
        signer.epid_sign(msg, &sig, 1);
        auto end_sign = std::chrono::high_resolution_clock::now();

        auto sign_duration =
            std::chrono::duration_cast<std::chrono::microseconds>(end_sign -
                                                                  start_sign);
        double sign_time_ms = sign_duration.count() / 1000.0;
        sign_times.push_back(sign_time_ms);

        auto start_verify = std::chrono::high_resolution_clock::now();
        bool verify_result = signer.epid_verify(msg, &sig);
        auto end_verify = std::chrono::high_resolution_clock::now();

        auto verify_duration =
            std::chrono::duration_cast<std::chrono::microseconds>(end_verify -
                                                                  start_verify);
        double verify_time_ms = verify_duration.count() / 1000.0;
        verify_times.push_back(verify_time_ms);

        std::cout << "  Sign time: " << sign_time_ms
                  << " ms, Verify time: " << verify_time_ms
                  << " ms, Result: " << (verify_result ? "pass" : "fail")
                  << std::endl;
    }

    double avg_sign_time =
        std::accumulate(sign_times.begin(), sign_times.end(), 0.0) /
        sign_times.size();
    double avg_verify_time =
        std::accumulate(verify_times.begin(), verify_times.end(), 0.0) /
        verify_times.size();

    std::cout << "\n=== Performance Test Results ===" << std::endl;
    std::cout << "Iterations: " << iterations << std::endl;
    std::cout << "Average Sign Time: " << avg_sign_time << " ms" << std::endl;
    std::cout << "Average Verify Time: " << avg_verify_time << " ms"
              << std::endl;
    std::cout << "Total Average Time per Test: "
              << (avg_sign_time + avg_verify_time) << " ms" << std::endl;
    std::cout << "=================================" << std::endl;
}
int main() {
    std::vector<uint8_t> msg = {0x11, 0x22};

    // Run single test case
    // std::cout << "=== Single Test Case ===" << std::endl;
    // test_case(msg, 128, 8, 8, 0, 16, 16, 1024*1024, 100);

    // Run performance test with 10 iterations
    std::cout << "\n=== Performance Test ===" << std::endl;
    performance_test(msg, 128, 8, 8, 0, 16, 16, 1024*1024, 100, 10);

    return 0;
}