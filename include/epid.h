#include <openssl/aes.h>
#include <openssl/evp.h>

#include <chrono>

#include "aes.h"
#include "aes_prove.h"
#include "hamming_weight_prove.h"
#include "instances.h"
#include "parameters.h"
#include "random_oracle.h"
#include "randomness.h"
#include "universal_hashing.h"
#include "utils.h"
#include "vole.h"
struct signature_t {
    std::vector<uint8_t> iv;
    std::vector<uint8_t> c;
    std::vector<uint8_t> u_tilde;
    std::vector<uint8_t> A_1_tilde_bytes;
    std::vector<uint8_t> A_2_tilde_bytes;
    std::vector<uint8_t> d;
    std::vector<uint8_t> chall_3;
    std::vector<uint8_t*> pdec;
    std::vector<uint8_t*> com;

    std::vector<uint8_t> t;
    std::vector<uint8_t> r;

    // // for debug
    // std::vector<uint8_t> chall_2;
    // std::vector<uint8_t> chall_1;
    // std::vector<uint8_t> A_0_tilde_bytes;
    // std::vector<bf128_t> a_0_vec;
    // std::vector<bf128_t> a_1_vec;
    // std::vector<bf128_t> a_2_vec;
    // std::vector<bf128_t> a_0_vec_test;
    // std::vector<bf128_t> a_1_vec_test;
    // std::vector<uint8_t> witness;
    

    // Destructor to clean up dynamically allocated memory
    ~signature_t() {
        for (auto* ptr : pdec) {
            delete[] ptr;
        }
        for (auto* ptr : com) {
            delete[] ptr;
        }
    }
};

class epid {  // a class to simulate the sign / verify process of our EPID
              // scheme
   public:
    int lambda_;
    int lambda_bytes_;
    const int iv_size_ = 16;

   public:
    epid(int lambda, int k0, int k1, int tau0, int tau1, int tau,
         unsigned int member_num, unsigned int test_srl_size)
        : lambda_(lambda), lambda_bytes_(lambda / 8) {
        params_.lambda = lambda;
        params_.k1 = k1;
        params_.k0 = k0;
        params_.tau0 = tau0;
        params_.tau1 = tau1;
        params_.tau = tau;
        member_num_ = member_num;
        // init group
        auto join_start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < member_num_; i++) {
            issue_join();
        }
        auto join_end = std::chrono::high_resolution_clock::now();
        std::cout << "issue_join loop time: "
                  << std::chrono::duration<double, std::milli>(join_end -
                                                               join_start)
                         .count()
                  << " ms" << std::endl;

        auto tree_start = std::chrono::high_resolution_clock::now();
        gen_tree();
        auto tree_end = std::chrono::high_resolution_clock::now();
        std::cout << "gen_tree time: "
                  << std::chrono::duration<double, std::milli>(tree_end -
                                                               tree_start)
                         .count()
                  << " ms" << std::endl;

        // Generate SRL with some revoked signatures
        auto srl_start = std::chrono::high_resolution_clock::now();
        generate_srl(test_srl_size);  // Generate test revoked signatures
        auto srl_end = std::chrono::high_resolution_clock::now();
        std::cout << "generate_srl time: "
                  << std::chrono::duration<double, std::milli>(srl_end -
                                                               srl_start)
                         .count()
                  << " ms" << std::endl;
    }

    void issue_join();
    void issue_rejoin();

    void epid_sign(const std::vector<uint8_t>& msg, signature_t* sig,
                   unsigned int signer_index);
    bool epid_verify(const std::vector<uint8_t>& msg, const signature_t* sig);

    void key_revoke();
    void sig_revoke();

    void generate_srl(unsigned int srl_size);

   private:
    void hash_mu(const std::vector<uint8_t>& msg, std::vector<uint8_t>& mu);

    void hash_challenge_1(const std::vector<uint8_t>& mu,
                          const std::vector<uint8_t>& hcom,
                          const std::vector<uint8_t>& c,
                          const std::vector<uint8_t>& iv,
                          std::vector<uint8_t>& chall_1, unsigned int ell,
                          unsigned int tau);
    void hash_challenge_2(std::vector<uint8_t>& chall_2,
                          const std::vector<uint8_t>& chall_1,
                          const std::vector<uint8_t>& u_tilde,
                          const std::vector<uint8_t>& h_v,
                          const std::vector<uint8_t>& d, unsigned int lambda,
                          unsigned int ell);
    void hash_challenge_3(std::vector<uint8_t>& chall_3,
                          const std::vector<uint8_t>& chall_2,
                          const std::vector<uint8_t>& a_0_tilde,
                          const std::vector<uint8_t>& a_1_tilde,
                          const std::vector<uint8_t>& a_2_tilde,
                          unsigned int lambda);

    // GF2_128 zk_hash(const std::vector<uint8_t>& sd,
    //                 const std::vector<GF2_128>& x_0, GF2_128& x_1);

    void gen_rootkey_iv(const std::vector<uint8_t>& mu, unsigned int index,
                        std::vector<uint8_t>& rootkey,
                        std::vector<uint8_t>& iv);

    void hash_1_witness_multiplexer(const std::vector<uint8_t>& input_0,
                                    const std::vector<uint8_t>& input_1,
                                    const std::vector<uint8_t>& s_byte,
                                    uint8_t* witness);

    void gen_witness_merkle_tree(uint8_t* witness, unsigned int signer_index);

    bf128_t zk_hash(const std::vector<uint8_t>& sd,
                    const std::vector<bf128_t>& x_0, const bf128_t& x_1);

    void gen_sk();
    void gen_challenge();
    void gen_x();
    void cal_t(const std::vector<uint8_t>& sk, const std::vector<uint8_t>& c_i);
    void gen_tree();

   private:
    int ell_ = 0;
    int member_num_ = 0;
    paramset_t params_;

    std::vector<std::vector<uint8_t>> tree_;
    std::vector<std::vector<uint8_t>> skey_set;
    std::vector<std::vector<uint8_t>> challenge_set;
    std::vector<std::vector<uint8_t>> t_set;
    std::vector<std::vector<uint8_t>> x_set;

    // SRL (Signature Revocation List)
    std::vector<signature_t> srl;

    const std::vector<uint8_t> s_0_ = {0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                                       0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                                       0x00, 0x00, 0x00, 0x00};
    const std::vector<uint8_t> s_1_ = {0x01, 0x00, 0x00, 0x00, 0x00, 0x00,
                                       0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                                       0x00, 0x00, 0x00, 0x00};
};
