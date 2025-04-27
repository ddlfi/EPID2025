#include <chrono>
#include <random>

#include "commit_prove.h"
#include "hamming_weight_prove.h"
#include "product_proof.h"
#include "signature.h"

struct baps_signature_t {
    signature_t* sig_attribute;
    signature_t* sig_policy;
    std::vector<uint8_t> iv;
    std::vector<uint8_t> c;
    std::vector<uint8_t> u_tilde;
    std::vector<uint8_t> A_1_tilde_bytes;
    std::vector<uint8_t> d;
    std::vector<uint8_t> chall_3;
    std::vector<uint8_t*> pdec;
    std::vector<uint8_t*> com;

    /////for debug

    // std::vector<u_int8_t> chall_2;
    // std::vector<u_int8_t> chall_1;
    // std::vector<field::GF2_128> A_0;
    // std::vector<field::GF2_128> A_1;

    baps_signature_t() {
        sig_attribute = new signature_t();
        sig_policy = new signature_t();
    }

    ~baps_signature_t() {
        delete sig_attribute;
        delete sig_policy;
    }

    baps_signature_t(const baps_signature_t&) = delete;
    baps_signature_t& operator=(const baps_signature_t&) = delete;
};

struct baps_attestation_t {
    std::vector<uint8_t> iv;
    std::vector<uint8_t> c;
    std::vector<uint8_t> u_tilde;
    std::vector<uint8_t> A_1_tilde_bytes;
    std::vector<uint8_t> d;
    std::vector<uint8_t> chall_3;
    std::vector<uint8_t*> pdec;
    std::vector<uint8_t*> com;

    /////for debug

    // std::vector<u_int8_t> chall_2;
    // std::vector<u_int8_t> chall_1;
    // std::vector<field::GF2_128> A_0;
    // std::vector<field::GF2_128> A_1;
    // std::vector<field::GF2_128> v;
};

class Baps {
    Signature attribute;
    Signature policy;

    const int lambda_ = 128;
    const int lambda_bytes_ = 16;
    const int iv_size_ = 16;
    const int k_ = 384;

   private:
    void hash_mu_baps(std::vector<uint8_t>& mu,
                      const std::vector<uint8_t>& msg);
    void hash_mu_disclose_baps(std::vector<uint8_t>& mu);
    void gen_rootkey_iv_baps(const std::vector<uint8_t>& mu,
                             std::vector<uint8_t>& rootkey,
                             std::vector<uint8_t>& iv);

    void get_y_baps();

    void get_z_baps();

    void get_x_z_baps();

    void get_b_baps();

    void gen_G_baps();

    void gen_w_baps();

    void calculate_t_baps();

   public:
    Baps(int num_)
        : attribute(num_),
          policy(2 * num_),
          G1_(k_, std::vector<uint8_t>(k_ * k_)),
          G2_(k_, std::vector<uint8_t>(k_)),
          t_(k_) {
        num_members_ = num_;
        attribute.get_x_(attribute_, attribute_index, 1);
        policy.get_x_(policy_, policy_index, 2);
        gen_w_baps();
        get_y_baps();
        get_z_baps();
        get_x_z_baps();  // x||z
        get_b_baps();
        gen_G_baps();
        calculate_t_baps();
        params_.lambda = 128;
        params_.k1 = 8;
        params_.k0 = 8;
        params_.tau0 = 0;
        params_.tau1 = 16;
        params_.tau = 16;
    }

    void baps_sign(int index_attribute, int index_policy,
                   const std::vector<uint8_t>& msg, baps_signature_t* sig,double& during_time);
    bool baps_verify(const std::vector<uint8_t>& msg, baps_signature_t* sig,double& during_time);

    void baps_disclose(baps_attestation_t* attestation,double& during_time);

    bool baps_disclose_verify(baps_attestation_t* attestation,double& during_time);

    void get_constrain_prover_disclose(const field::GF2_128* v,
                                       const uint8_t* witness,
                                       field::GF2_128* A_0,
                                       field::GF2_128* A_1);

    void get_constrain_verifier_disclose(const field::GF2_128* q,
                                         field::GF2_128 delta,
                                         field::GF2_128* B);

    void baps_commit(const std::vector<uint8_t>& msg,
                     const std::vector<uint8_t>& r, uint8_t* witness);

   private:
    const int attribute_index = 0;
    const int policy_index = 0;
    int num_members_;

    std::vector<uint8_t> attribute_;
    std::vector<uint8_t> policy_;
    std::vector<uint8_t> w_;
    std::vector<uint8_t> y_;
    std::vector<uint8_t> z_;
    std::vector<uint8_t> x_z;  // x||z
    std::vector<uint8_t> b_;   // vp||x
    std::vector<uint8_t> t_;

    std::vector<std::vector<uint8_t>> G1_;
    std::vector<std::vector<uint8_t>> G2_;
    paramset_t params_;
};