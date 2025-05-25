#include "hamming_weight_prove.h"
#include "linear_prove.h"
#include "random_oracle.h"
#include "randomness.h"
#include "universal_hashing.h"
#include "utils.h"
#include "vole.h"

template <typename GF2>
struct signature_t {
    std::vector<uint8_t> iv;
    std::vector<uint8_t> c;
    std::vector<uint8_t> u_tilde;
    std::vector<uint8_t> A_1_tilde_bytes;
    std::vector<uint8_t> d;
    std::vector<uint8_t> chall_3;
    std::vector<uint8_t*> pdec;
    std::vector<uint8_t*> com;

    // for verify hamming weight
    std::vector<uint8_t> hamming_weight_vec;
    std::vector<GF2> hamming_weight_v;

    // for verify linear operation:
    std::vector<GF2> A_linear;
};

template <typename GF2>
class Signature {
   public:
    int lambda_;
    int lambda_bytes_;
    const int iv_size_ = 16;
    int m_;
    int k_;  // 矩阵行数
    int hamming_weight_;

   public:
    Signature(int lambda, int k0, int k1, int tau0, int tau1, int tau, int m,
              int k, int hamming_weight)
        : matrix_(k_, std::vector<uint8_t>((m_ + 7) / 8)),
          x_((m_ + 7) / 8),
          lambda_(lambda),
          lambda_bytes_(lambda / 8),
          m_(m),
          k_(k),
          hamming_weight_(hamming_weight) {
        params_.lambda = lambda;
        params_.k1 = k1;
        params_.k0 = k0;
        params_.tau0 = tau0;
        params_.tau1 = tau1;
        params_.tau = tau;

        calculate_params();
        gen_matrix();
        gen_x();
        calculate_y();
    }
    void sign(const std::vector<uint8_t>& msg, signature_t<GF2>* sig);
    bool verify(const std::vector<uint8_t>& msg, const signature_t<GF2>* sig);

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

    GF2 zk_hash(const std::vector<uint8_t>& sd, const std::vector<GF2>& x_0,
                GF2& x_1);

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
