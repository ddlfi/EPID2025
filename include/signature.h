#include "hamming_weight_prove.h"
#include "random_oracle.h"
#include "randomness.h"
#include "universal_hashing.h"
#include "utils.h"
#include "vole.h"
#include "linear_prove.h"

struct signature_t {
    std::vector<uint8_t> iv;
    std::vector<uint8_t> c;
    std::vector<uint8_t> u_tilde;
    std::vector<uint8_t> A_1_tilde_bytes;
    std::vector<uint8_t> d;
    std::vector<uint8_t> chall_3;
    std::vector<uint8_t*> pdec;
    std::vector<uint8_t*> com;

    //for verify hamming weight
    std::vector<uint8_t> hamming_weight_vec;
    std::vector<field::GF2_128> hamming_weight_v;

    //for verify linear operation:
    std::vector<field::GF2_128> A_linear;
    
    // /////for debug

    // std::vector<u_int8_t> chall_2;
    // std::vector<field::GF2_128> A_0;
    // std::vector<field::GF2_128> A_1;

    // std::vector<uint8_t> witness;
    // std::vector<field::GF2_128> v_combined;

    // std::vector<field::GF2_128> v_field;
};

void hash_challenge_1(const std::vector<uint8_t>& mu,
                      const std::vector<uint8_t>& hcom,
                      const std::vector<uint8_t>& c,
                      const std::vector<uint8_t>& iv,
                      std::vector<uint8_t>& chall_1, unsigned int lambda,
                      unsigned int ell, unsigned int tau);
void hash_challenge_2(std::vector<uint8_t>& chall_2,
                      const std::vector<uint8_t>& chall_1,
                      const std::vector<uint8_t>& u_tilde,
                      const std::vector<uint8_t>& h_v,
                      const std::vector<uint8_t>& d, unsigned int lambda,
                      unsigned int ell);
void hash_challenge_3(std::vector<uint8_t>& chall_3,
                      const std::vector<uint8_t>& chall_2,
                      const std::vector<uint8_t>& a_tilde,
                      const std::vector<uint8_t>& b_tilde, unsigned int lambda);

field::GF2_128 zk_hash(const std::vector<uint8_t>& sd,
                       const std::vector<field::GF2_128>& x_0,
                       field::GF2_128& x_1);

class Signature {
   public:
    const int lambda_ = 128;
    const int lambda_bytes_ = lambda_ / 8;
    const int iv_size_ = 16;
    const int m_ = 1268;
    const int k_ = 634;  // 矩阵行数
    const int hamming_weight = 133;

   public:
    Signature()
        : matrix_(k_, std::vector<uint8_t>((m_ + 7) / 8)), x_((m_ + 7) / 8) {
        params_.lambda = 128;
        params_.k1 = 8;
        params_.k0 = 8;
        params_.tau0 = 0;
        params_.tau1 = 16;
        params_.tau = 16;

        calculate_params();
        gen_matrix();
        gen_x();
        calculate_y();
    }
    void sign(const std::vector<uint8_t>& msg, signature_t* sig);
    bool verify(const std::vector<uint8_t>& msg, const signature_t* sig);

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
                          const std::vector<uint8_t>& a_tilde,
                          const std::vector<uint8_t>& b_tilde,
                          unsigned int lambda);

    field::GF2_128 zk_hash(const std::vector<uint8_t>& sd,
                           const std::vector<field::GF2_128>& x_0,
                           field::GF2_128& x_1);

    void gen_rootkey_iv(const std::vector<uint8_t>& mu,
                        std::vector<uint8_t>& rootkey,
                        std::vector<uint8_t>& iv);

    void calculate_params();

    void gen_matrix();
    void gen_x();
    void calculate_y();

   private:
    int ell_ = 0;
    int muti_num_ = 0;
    std::vector<std::vector<uint8_t>> matrix_;
    std::vector<uint8_t> x_;
    std::vector<uint8_t> y_;
    paramset_t params_;
};
