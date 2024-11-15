#include "path_prove.h"
#include "rain.h"
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
    std::vector<uint8_t> d;
    std::vector<uint8_t> chall_3;
    std::vector<uint8_t*> pdec;
    std::vector<uint8_t*> com;

    // /////for debug

    // std::vector<u_int8_t> chall_2;
    // std::vector<field::GF2_128> A_0;
    // std::vector<field::GF2_128> A_1;

    // std::vector<uint8_t> witness;
    // std::vector<field::GF2_128> v_combined;

    // std::vector<field::GF2_128> v_field;
};

class Signature {
   public:
    int key_num_;
    int tree_node_num_;
    const int lambda_ = 128;
    const int lambda_bytes_ = 16;
    const int iv_size_ = 16;

   public:
    Signature(unsigned int key_num)
        : skey_(key_num, std::vector<uint8_t>(lambda_bytes_)),
          pkey_(key_num, std::vector<uint8_t>(lambda_bytes_)),
          tree_(2 * key_num, std::vector<uint8_t>(lambda_bytes_)) {
        key_num_ = key_num;
        tree_node_num_ = 2 * key_num;
        gen_skey();
        gen_pkey();
        gen_tree();
        params_.lambda = 128;
        params_.k1 = 8;
        params_.k0 = 8;
        params_.tau0 = 0;
        params_.tau1 = 16;
        params_.tau = 16;
    }
    void sign(unsigned int signer_index, const std::vector<uint8_t>& msg,
              signature_t* sig);
    bool verify(const std::vector<uint8_t>& msg, const signature_t* sig);

   private:
    void gen_pkey();

    void gen_tree();

    void gen_skey();

    void hash_1(const std::vector<uint8_t>& input_0,
                const std::vector<uint8_t>& input_1,
                const std::vector<uint8_t>& s_byte,
                std::vector<uint8_t>& output, uint8_t* witness, bool flag);

    void hash_pk_msg(const std::vector<uint8_t>& msg, std::vector<uint8_t>& mu);

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
                        const uint8_t signer_index,
                        std::vector<uint8_t>& rootkey,
                        std::vector<uint8_t>& iv);

    void gen_witness(uint8_t* witness, unsigned int index);

   private:
    const std::vector<uint8_t> rain_msg_ = {0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                                            0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                                            0x00, 0x00, 0x00, 0x00};

    const std::vector<uint8_t> s_0_ = {0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                                       0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                                       0x00, 0x00, 0x00, 0x00};
    const std::vector<uint8_t> s_1_ = {0x01, 0x00, 0x00, 0x00, 0x00, 0x00,
                                       0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                                       0x00, 0x00, 0x00, 0x00, 0x00};

    std::vector<std::vector<uint8_t>> skey_;

    std::vector<std::vector<uint8_t>> pkey_;

    std::vector<std::vector<uint8_t>> tree_;

    paramset_t params_;
};
