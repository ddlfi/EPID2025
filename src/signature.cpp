#include "signature.h"

#include <chrono>
#include <cmath>
#include <random>

void Signature::gen_skey() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint8_t> dis(0, 255);

    for (auto& skey_i : skey_) {
        for (auto& value : skey_i) {
            value = dis(gen);
        }
    }
}

void Signature::gen_pkey() {
    for (auto i = 0; i < key_num_; i++) {
        rain(skey_[i], rain_msg_, pkey_[i], nullptr, 0);
    }
}

void Signature::gen_tree() {
    for (int i = key_num_; i < tree_node_num_; i++) {
        tree_[i] = pkey_[i - key_num_];
    }
    int index = key_num_ - 1;
    while (index != 0) {
        hash_1(tree_[2 * index], tree_[2 * index + 1], s_0_, tree_[index],
               nullptr, 0);
        index--;
    }
}

void Signature::hash_1(const std::vector<uint8_t>& input_0,
                       const std::vector<uint8_t>& input_1,
                       const std::vector<uint8_t>& s_byte,
                       std::vector<uint8_t>& output, uint8_t* witness,
                       bool flag) {
    if (flag && witness) {
        memcpy(witness, input_0.data(), lambda_bytes_);
        memcpy(witness + lambda_bytes_, input_1.data(), lambda_bytes_);
        memcpy(witness + lambda_bytes_ * 2, s_byte.data(), lambda_bytes_);
        witness += lambda_bytes_ * 3;
    }
    field::GF2_128 i_0, i_1, s, s_, key, msg, rain_output, result;

    i_0.from_bytes(input_0.data());
    i_1.from_bytes(input_1.data());
    s.from_bytes(s_byte.data());
    s_ = (i_0 + i_1) * s;
    if (flag && witness) {
        s_.to_bytes(witness);
        witness += lambda_bytes_;
    }
    key = s_ + i_1;
    msg = s_ + i_0;

    rain(key, msg, rain_output, witness, flag);
    witness += lambda_bytes_ * 2;

    result = rain_output + msg;
    result.to_bytes(output.data());

    if (flag && witness) {
        result.to_bytes(witness);
    }
}

void Signature::hash_pk_msg(const std::vector<uint8_t>& msg,
                            std::vector<uint8_t>& mu) {
    H1_context_t h1_ctx;
    H1_init(&h1_ctx, lambda_);
    for (auto& pk : pkey_) {
        H1_update(&h1_ctx, pk.data(), pk.size());
    }
    H1_update(&h1_ctx, msg.data(), msg.size());
    H1_final(&h1_ctx, mu.data(), 2 * lambda_bytes_);
}

void Signature::hash_challenge_1(const std::vector<uint8_t>& mu,
                                 const std::vector<uint8_t>& hcom,
                                 const std::vector<uint8_t>& c,
                                 const std::vector<uint8_t>& iv,
                                 std::vector<uint8_t>& chall_1,
                                 unsigned int ell, unsigned int tau) {
    const unsigned int ell_hat_bytes =
        ell / 8 + lambda_bytes_ * 2 + UNIVERSAL_HASH_B;
    H2_context_t h2_ctx;
    H2_init(&h2_ctx, lambda_);
    H2_update(&h2_ctx, mu.data(), lambda_bytes_ * 2);
    H2_update(&h2_ctx, hcom.data(), lambda_bytes_ * 2);
    H2_update(&h2_ctx, c.data(), ell_hat_bytes * (tau - 1));
    H2_update(&h2_ctx, iv.data(), IV_SIZE);
    H2_final(&h2_ctx, chall_1.data(), 5 * lambda_bytes_ + 8);
}

void Signature::hash_challenge_2(std::vector<uint8_t>& chall_2,
                                 const std::vector<uint8_t>& chall_1,
                                 const std::vector<uint8_t>& u_tilde,
                                 const std::vector<uint8_t>& h_v,
                                 const std::vector<uint8_t>& d,
                                 unsigned int lambda, unsigned int ell) {
    const unsigned int lambda_bytes = lambda / 8;
    const unsigned int ell_bytes = ell / 8;
    const unsigned int u_tilde_bytes = lambda_bytes + UNIVERSAL_HASH_B;

    H2_context_t h2_ctx_1;
    H2_init(&h2_ctx_1, lambda);
    H2_update(&h2_ctx_1, chall_1.data(), 5 * lambda_bytes + 8);
    H2_update(&h2_ctx_1, u_tilde.data(), u_tilde_bytes);
    H2_update(&h2_ctx_1, h_v.data(), 2 * lambda_bytes);
    H2_update(&h2_ctx_1, d.data(), ell_bytes);
    H2_final(&h2_ctx_1, chall_2.data(), 3 * lambda_bytes + 8);
}

void Signature::hash_challenge_3(std::vector<uint8_t>& chall_3,
                                 const std::vector<uint8_t>& chall_2,
                                 const std::vector<uint8_t>& a_tilde,
                                 const std::vector<uint8_t>& b_tilde,
                                 unsigned int lambda) {
    const unsigned int lambda_bytes = lambda / 8;

    H2_context_t h2_ctx_2;
    H2_init(&h2_ctx_2, lambda);
    H2_update(&h2_ctx_2, chall_2.data(), 3 * lambda_bytes + 8);
    H2_update(&h2_ctx_2, a_tilde.data(), lambda_bytes);
    H2_update(&h2_ctx_2, b_tilde.data(), lambda_bytes);
    H2_final(&h2_ctx_2, chall_3.data(), lambda_bytes);
}

field::GF2_128 Signature::zk_hash(const std::vector<uint8_t>& sd,
                                  const std::vector<field::GF2_128>& x_0,
                                  field::GF2_128& x_1) {
    field::GF2_128 r_0, r_1, s, h_0, h_1;
    r_0.from_bytes(sd.data());
    r_1.from_bytes(sd.data() + lambda_bytes_);
    s.from_bytes(sd.data() + lambda_bytes_ * 2);

    uint64_t tmp;
    memcpy(&tmp, sd.data() + lambda_bytes_ * 3, 8UL);
    field::GF2_128 t(tmp);
    field::GF2_128 s_muti = s;
    field::GF2_128 t_muti = t;
    for (auto& x_0_i : x_0) {
        h_0 += s_muti * x_0_i;
        h_1 += t_muti * x_0_i;
        s_muti *= s;
        t_muti *= t;
    }
    return r_0 * h_0 + r_1 * h_1 + x_1;
}

void Signature::gen_rootkey_iv(const std::vector<uint8_t>& mu,
                               const uint8_t signer_index,
                               std::vector<uint8_t>& rootkey,
                               std::vector<uint8_t>& iv) {
    H3_context_t h3_ctx;
    H3_init(&h3_ctx, lambda_);
    H3_update(&h3_ctx, skey_[signer_index].data(), lambda_bytes_);
    H3_update(&h3_ctx, mu.data(), lambda_bytes_ * 2);
    H3_final(&h3_ctx, rootkey.data(), lambda_bytes_, iv.data());
}

void Signature::gen_witness(uint8_t* witness, unsigned int index) {
    std::vector<uint8_t> tmp;
    rain(skey_[index], rain_msg_, tmp, witness, 1);
    witness +=
        3 * lambda_bytes_;  // 最后一个witness和下一次hash重叠了，因此+3而非4
    index = index + key_num_;
    while (index != 1) {
        if (index % 2) {
            hash_1(tree_[index], tree_[index - 1], s_1_, tmp, witness, 1);
            witness += 6 * lambda_bytes_;
            index = index / 2;
        } else {
            hash_1(tree_[index], tree_[index + 1], s_0_, tmp, witness, 1);
            witness += 6 * lambda_bytes_;
            index = index / 2;
        }
    }
}

void Signature::sign(unsigned int signer_index, const std::vector<uint8_t>& msg,
                     signature_t* sig) {
    const unsigned int ell = (4 + 6 * log2(key_num_)) * lambda_;
    const unsigned int muti_times = 3 + 4 * log2(key_num_);
    const unsigned int ell_bytes = ell / 8;
    const unsigned int ell_hat = ell + lambda_ * 2 + UNIVERSAL_HASH_B_BITS;
    const unsigned int ell_hat_bytes = ell_hat / 8;

    // auto start_time = std::chrono::high_resolution_clock::now();

    std::vector<uint8_t> mu(2 * lambda_bytes_);
    hash_pk_msg(msg, mu);

    std::vector<uint8_t> rootkey(lambda_bytes_);
    sig->iv.resize(iv_size_);
    // std::vector<uint8_t> iv(iv_size_);
    gen_rootkey_iv(mu, signer_index, rootkey, sig->iv);

    auto vole_commit_start_time = std::chrono::high_resolution_clock::now();

    std::vector<uint8_t> hcom(lambda_bytes_ * 2);
    std::vector<vec_com_t> vecCom(params_.tau);
    std::vector<uint8_t> u(ell_hat_bytes);
    std::vector<uint8_t*> V(lambda_);
    V[0] = new uint8_t[lambda_ * ell_hat_bytes];
    for (unsigned int i = 1; i < lambda_; ++i) {
        V[i] = V[0] + i * ell_hat_bytes;
    }

    sig->c.resize((params_.tau - 1) * ell_hat_bytes);
    // std::vector<uint8_t> c((params_.tau - 1) * ell_hat_bytes);

    vole_commit(rootkey.data(), sig->iv.data(), ell_hat, &(params_),
                hcom.data(), vecCom.data(), sig->c.data(), u.data(), V.data());

    auto vole_commit_end_time = std::chrono::high_resolution_clock::now();

    std::vector<uint8_t> chall_1(5 * lambda_bytes_ + 8);
    hash_challenge_1(mu, hcom, sig->c, sig->iv, chall_1, ell, params_.tau);

    sig->u_tilde.resize(lambda_bytes_ + UNIVERSAL_HASH_B);
    // std::vector<uint8_t> u_tilde(lambda_bytes_ + UNIVERSAL_HASH_B);
    vole_hash(sig->u_tilde.data(), chall_1.data(), u.data(), ell, lambda_);

    std::vector<uint8_t> h_v(lambda_bytes_ * 2);
    {
        H1_context_t h1_ctx_1;
        H1_init(&h1_ctx_1, lambda_);
        std::vector<uint8_t> V_tilde(lambda_bytes_ + UNIVERSAL_HASH_B);
        for (unsigned int i = 0; i != lambda_; ++i) {
            // Step 7
            vole_hash(V_tilde.data(), chall_1.data(), V[i], ell, lambda_);
            // Step 8
            H1_update(&h1_ctx_1, V_tilde.data(),
                      lambda_bytes_ + UNIVERSAL_HASH_B);
        }
        // Step: 8
        H1_final(&h1_ctx_1, h_v.data(), lambda_bytes_ * 2);
    }

    std::vector<uint8_t> witness(ell_bytes);
    gen_witness(witness.data(), signer_index);
    sig->d.resize(ell_bytes);
    // std::vector<uint8_t> d(ell_bytes);
    xor_u8_array(witness.data(), u.data(), sig->d.data(), ell_bytes);

    std::vector<uint8_t> chall_2(3 * lambda_bytes_ + 8);
    hash_challenge_2(chall_2, chall_1, sig->u_tilde, h_v, sig->d, lambda_, ell);

    std::vector<field::GF2_128> v_gf_128_vec(ell_hat);
    std::vector<field::GF2_128> v_combined_gf_128_vec(ell_hat / lambda_);
    convert_vec_to_field(V.data(), v_gf_128_vec.data(), ell_hat, lambda_);
    gen_combined_field_vec(v_gf_128_vec.data(), v_combined_gf_128_vec.data(),
                           ell_hat);

    std::vector<field::GF2_128> A_0(muti_times);
    std::vector<field::GF2_128> A_1(muti_times);

    path_prove(witness.data(), v_combined_gf_128_vec.data(),
               v_gf_128_vec.data(), rain_msg_, A_0.data(), A_1.data(),
               log2(key_num_));

    field::GF2_128 u_star;
    u_star.from_bytes(u.data() + ell_bytes);
    field::GF2_128 A_0_tilde =
        zk_hash(chall_2, A_0, v_combined_gf_128_vec[ell / lambda_]);
    field::GF2_128 A_1_tilde = zk_hash(chall_2, A_1, u_star);

    std::vector<uint8_t> A_0_tilde_bytes(lambda_bytes_);
    sig->A_1_tilde_bytes.resize(lambda_bytes_);
    A_0_tilde.to_bytes(A_0_tilde_bytes.data());
    A_1_tilde.to_bytes(sig->A_1_tilde_bytes.data());

    sig->chall_3.resize(lambda_bytes_);
    // std::vector<uint8_t> chall_3(lambda_bytes_);
    sig->pdec.resize(params_.tau);
    sig->com.resize(params_.tau);
    // std::vector<uint8_t*> pdec(params_.tau);
    // std::vector<uint8_t*> com(params_.tau);

    hash_challenge_3(sig->chall_3, chall_2, sig->A_1_tilde_bytes,
                     A_0_tilde_bytes, lambda_);

    for (unsigned int i = 0; i < params_.tau; i++) {
        // Step 20
        uint8_t s_[12];
        ChalDec(sig->chall_3.data(), i, params_.k0, params_.tau0, params_.k1,
                params_.tau1, s_);

        // Step 21
        const unsigned int depth = i < params_.tau0 ? params_.k0 : params_.k1;
        sig->pdec[i] = new uint8_t[depth * lambda_bytes_];
        sig->com[i] = new uint8_t[2 * lambda_bytes_];
        vector_open(vecCom[i].k, vecCom[i].com, s_, sig->pdec[i], sig->com[i],
                    depth, lambda_bytes_);
        vec_com_clear(&vecCom[i]);
    }

    auto end_time = std::chrono::high_resolution_clock::now();

    // auto total_time = end_time - start_time;
    auto vole_time = vole_commit_end_time - vole_commit_start_time;

    // std::cout << "sign total time is : "
    //           << std::chrono::duration<double,
    //           std::milli>(total_time).count()
    //           << " ms" << std::endl;
    std::cout << "vole commit time is : "
              << std::chrono::duration<double, std::milli>(vole_time).count()
              << " ms" << std::endl;
    delete V[0];
}

bool Signature::verify(const std::vector<uint8_t>& msg,
                       const signature_t* sig) {
    const unsigned int ell = (4 + 6 * log2(key_num_)) * lambda_;
    const unsigned int muti_times = 3 + 4 * log2(key_num_);
    const unsigned int ell_bytes = ell / 8;
    const unsigned int ell_hat = ell + lambda_ * 2 + UNIVERSAL_HASH_B_BITS;
    const unsigned int ell_hat_bytes = ell_hat / 8;

    // auto start_time = std::chrono::high_resolution_clock::now();

    std::vector<uint8_t> mu(2 * lambda_bytes_);
    hash_pk_msg(msg, mu);

    auto vole_reconstruct_start_time =
        std::chrono::high_resolution_clock::now();

    std::vector<uint8_t*> Q(lambda_);
    std::vector<uint8_t> hcom(lambda_bytes_ * 2);
    Q[0] = new uint8_t[lambda_ * ell_hat_bytes];
    for (unsigned int i = 1; i < lambda_; ++i) {
        Q[i] = Q[0] + i * ell_hat_bytes;
    }
    vole_reconstruct(sig->iv.data(), sig->chall_3.data(), sig->pdec.data(),
                     sig->com.data(), hcom.data(), Q.data(), ell_hat, &params_);

    auto vole_reconstruct_end_time = std::chrono::high_resolution_clock::now();

    std::vector<uint8_t> chall_1(5 * lambda_bytes_ + 8);
    hash_challenge_1(mu, hcom, sig->c, sig->iv, chall_1, ell, params_.tau);

    std::vector<uint8_t*> Q_(lambda_);
    Q_[0] = new uint8_t[lambda_ * ell_hat_bytes];
    for (unsigned int i = 1; i < lambda_; ++i) {
        Q_[i] = Q_[0] + i * ell_hat_bytes;
    }

    std::vector<uint8_t*> Dtilde(lambda_);
    Dtilde[0] = new uint8_t[lambda_ * (lambda_bytes_ + UNIVERSAL_HASH_B)];
    for (unsigned int i = 1; i < lambda_; ++i) {
        Dtilde[i] = Dtilde[0] + i * (lambda_bytes_ + UNIVERSAL_HASH_B);
    }
    memset(Dtilde[0], 0, lambda_ * (lambda_bytes_ + UNIVERSAL_HASH_B));

    unsigned int Dtilde_idx = 0;
    unsigned int q_idx = 0;
    for (unsigned int i = 0; i < params_.tau; i++) {
        const unsigned int depth = i < params_.tau0 ? params_.k0 : params_.k1;

        // Step 11
        uint8_t delta[8];
        ChalDec(sig->chall_3.data(), i, params_.k0, params_.tau0, params_.k1,
                params_.tau1, delta);
        // Step 16
        for (unsigned int j = 0; j != depth; ++j, ++Dtilde_idx) {
            // for scan-build
            assert(Dtilde_idx < lambda_);
            masked_xor_u8_array(Dtilde[Dtilde_idx], sig->u_tilde.data(),
                                Dtilde[Dtilde_idx], delta[j],
                                lambda_bytes_ + UNIVERSAL_HASH_B);
        }

        if (i == 0) {
            // Step 8
            memcpy(Q_[q_idx], Q[q_idx], ell_hat_bytes * depth);
            q_idx += depth;
        } else {
            // Step 14
            for (unsigned int d = 0; d < depth; ++d, ++q_idx) {
                masked_xor_u8_array(Q[q_idx],
                                    sig->c.data() + (i - 1) * ell_hat_bytes,
                                    Q_[q_idx], delta[d], ell_hat_bytes);
            }
        }
    }

    std::vector<uint8_t> h_v(lambda_bytes_ * 2);
    {
        H1_context_t h1_ctx_1;
        H1_init(&h1_ctx_1, lambda_);
        std::vector<uint8_t> Q_tilde(lambda_bytes_ + UNIVERSAL_HASH_B);
        for (unsigned int i = 0; i < lambda_; i++) {
            vole_hash(Q_tilde.data(), chall_1.data(), Q_[i], ell, lambda_);
            xor_u8_array(Q_tilde.data(), Dtilde[i], Q_tilde.data(),
                         lambda_bytes_ + UNIVERSAL_HASH_B);
            H1_update(&h1_ctx_1, Q_tilde.data(),
                      lambda_bytes_ + UNIVERSAL_HASH_B);
        }
        H1_final(&h1_ctx_1, h_v.data(), lambda_bytes_ * 2);
    }
    delete Dtilde[0];

    std::vector<uint8_t> chall_2(3 * lambda_bytes_ + 8);
    hash_challenge_2(chall_2, chall_1, sig->u_tilde, h_v, sig->d, lambda_, ell);

    for (unsigned int i = 0, col = 0; i < params_.tau; i++) {
        unsigned int depth = i < params_.tau0 ? params_.k0 : params_.k1;
        uint8_t decoded_challenge[MAX_DEPTH];
        ChalDec(sig->chall_3.data(), i, params_.k0, params_.tau0, params_.k1,
                params_.tau1, decoded_challenge);
        for (unsigned int j = 0; j < depth; j++, ++col) {
            if (decoded_challenge[j] == 1) {
                xor_u8_array(sig->d.data(), Q_[col], Q_[col], ell_bytes);
            }
        }
    }

    std::vector<field::GF2_128> q_gf_128_vec(ell_hat);
    std::vector<field::GF2_128> q_combined_gf_128_vec(ell_hat / lambda_);
    convert_vec_to_field(Q_.data(), q_gf_128_vec.data(), ell_hat, lambda_);
    gen_combined_field_vec(q_gf_128_vec.data(), q_combined_gf_128_vec.data(),
                           ell_hat);

    field::GF2_128 delta_field;
    delta_field.from_bytes(sig->chall_3.data());

    std::vector<field::GF2_128> B(muti_times);
    path_verify(q_combined_gf_128_vec.data(), q_gf_128_vec.data(), delta_field,
                rain_msg_, B.data(), log2(key_num_));

    field::GF2_128 zero;
    field::GF2_128 B_tilde =
        zk_hash(chall_2, B, q_combined_gf_128_vec[ell / lambda_]);

    field::GF2_128 A_0_tilde, A_1_tilde;
    std::vector<uint8_t> A_0_tilde_bytes(lambda_bytes_);

    A_1_tilde.from_bytes(sig->A_1_tilde_bytes.data());
    A_0_tilde = B_tilde - delta_field * A_1_tilde;
    A_0_tilde.to_bytes(A_0_tilde_bytes.data());

    std::vector<uint8_t> chall_3(lambda_bytes_);
    hash_challenge_3(chall_3, chall_2, sig->A_1_tilde_bytes, A_0_tilde_bytes,
                     lambda_);

    // auto end_time = std::chrono::high_resolution_clock::now();

    // auto total_time = end_time - start_time;
    auto vole_time = vole_reconstruct_end_time - vole_reconstruct_start_time;

    // std::cout << "verify total time is : "
    //           << std::chrono::duration<double,
    //           std::milli>(total_time).count()
    //           << " ms" << std::endl;
    std::cout << "vole reconstruct time is : "
              << std::chrono::duration<double, std::milli>(vole_time).count()
              << " ms" << std::endl;

    delete Q[0];
    delete Q_[0];

    return memcmp(chall_3.data(), sig->chall_3.data(), lambda_bytes_) == 0;
}
