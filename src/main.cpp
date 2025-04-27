// #include <chrono>
// #include <cstdint>
// #include <iostream>

// #include "signature.h"

// void test_case(unsigned int member_num, unsigned int signer_index,
//                const std::vector<uint8_t>& msg) {
//     Signature signer(member_num);
//     std::cout<<"Start test case : "<<member_num<<" members :"<<std::endl;
//     std::cout<<std::endl;
//     // std::cout << "         sign stage :" << std::endl;
//     auto time_1 = std::chrono::high_resolution_clock::now();
//     signer.gen_tree();
//     auto time_2 = std::chrono::high_resolution_clock::now();
//     auto tree_time = time_2 - time_1;
//     // std::cout << "merkle tree generation time is : "
//     //           << std::chrono::duration<double,
//     std::milli>(tree_time).count()
//     //           << " ms" << std::endl;

//     signature_t sig;
//     signer.sign(signer_index, msg, &sig);
//     auto time_3 = std::chrono::high_resolution_clock::now();

//     auto sign_time = time_3 - time_1;
//     std::cout << "total sign time is : "
//               << std::chrono::duration<double, std::milli>(sign_time).count()
//               << " ms" << std::endl;

//     // std::cout << "         verify stage :" << std::endl;
//     // std::cout << "merkle tree generation time is : "
//     //           << std::chrono::duration<double,
//     std::milli>(tree_time).count()
//     //           << " ms" << std::endl;
//     std::chrono::time_point<std::chrono::high_resolution_clock> time_4;
//     if (signer.verify(msg, &sig)) {
//         time_4 = std::chrono::high_resolution_clock::now();
//         std::cout << "total verify time is : "
//                   << std::chrono::duration<double, std::milli>(time_4 -
//                   time_3 + time_2 - time_1)
//                          .count()
//                   << " ms" << std::endl;
//         std::cout << "test case : " << member_num << " members pass!"
//                   << std::endl;
//     } else {
//         time_4 = std::chrono::high_resolution_clock::now();
//         std::cout << "total sign time is : "
//                   << std::chrono::duration<double, std::milli>(time_4 -
//                   time_3)
//                          .count()
//                   << " ms" << std::endl;
//         std::cout << "test case : " << member_num << " members fail!"
//                   << std::endl;
//     }
//     std::cout << std::endl;
// }

// int main() {
//     gen_field_base_128(field_base_128);
//     std::vector<uint8_t> msg = {0x11, 0x22};

//     test_case(1024U, 512U, msg);
//     test_case(32768,0,msg);
// }

#include <chrono>
#include <cstdint>
#include <iostream>

#include "baps.h"

void test_case(const std::vector<uint8_t>& msg, int num, double& sign_time,
               double& verify_time, double& disclose_sign_time,
               double& disclose_verify_time) {
    Baps baps_(num);
    baps_signature_t sig;
    baps_attestation_t attestation;
    baps_.baps_sign(0, 0, msg, &sig,sign_time);

    if (baps_.baps_verify(msg, &sig, verify_time)) {
        // std::cout << "test case : " << num << " members sign pass!"
        //           << std::endl;
    } else {
        // std::cout << "test case : " << num << " members sign fail!"
        //           << std::endl;
    }
    baps_.baps_disclose(&attestation,disclose_sign_time);
    if (baps_.baps_disclose_verify(&attestation,disclose_verify_time)) {
        // std::cout << "test case : " << num << " members disclose pass!"
        //           << std::endl;
    } else {
        // std::cout << "test case : " << num << " members disclose fail!"
        //           << std::endl;
    }
    // std::cout << std::endl;
}

void test_average_time(const std::vector<uint8_t>& msg, int num) {
    std::cout << "N1 and N2: 2^" << log2(num) << std::endl;
    double sign_time, verify_time, disclose_sign_time, disclose_verify_time;
    double t1 = 0 , t2 = 0, t3 = 0, t4 = 0;
    for (int i = 0; i < 10; i++) {
        test_case(msg, num,sign_time,verify_time,disclose_sign_time,disclose_verify_time);
        t1 += sign_time;
        t2 += verify_time;
        t3 += disclose_sign_time;
        t4 += disclose_verify_time;
    }
    std::cout<<"sign time is : "<< double(t1 / 10)<<"ms"<<std::endl;
    std::cout<<"verify time is : "<< double(t2 / 10)<<"ms"<<std::endl;
    std::cout<<"disclose time is : "<< double(t3 / 10)<<"ms"<<std::endl;
    std::cout<<"disclose verify time is : "<< double(t4 / 10)<<"ms"<<std::endl;
    std::cout << std::endl;
}
int main() {
    gen_field_base_128(field_base_128);
    std::vector<uint8_t> msg = {0x11, 0x22};
    std::cout << "128 bit security:" << std::endl;
    std::cout << std::endl;
    test_average_time(msg, 64UL);
    test_average_time(msg, 128UL);
    test_average_time(msg, 256UL);
    test_average_time(msg, 512UL);
    test_average_time(msg, 1024UL);
    test_average_time(msg, 2048UL);
    test_average_time(msg, 4096UL);
    test_average_time(msg, 8192UL);
    test_average_time(msg, 16384UL);
    test_average_time(msg, 32768UL);
}