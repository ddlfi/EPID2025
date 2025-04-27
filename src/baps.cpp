#include "baps.h"

bool compare_uint8_arrays(const std::vector<uint8_t>& arr1,
                          const std::vector<uint8_t>& arr2) {
    return (arr1.size() == arr2.size()) &&
           std::equal(arr1.begin(), arr1.end(), arr2.begin());
}

void Baps::get_y_baps() {
    int index = 0;
    for (int i = 0; i < lambda_ - 1; i++) {
        if (GET_BIT(w_[i / 8], i % 8)) index++;
    }
    for (int i = 0; i < lambda_; i++) {
        if (GET_BIT(attribute_[i / 8], i % 8)) index++;
    }
    y_.resize(2 * lambda_bytes_);
    setbit(&y_[index / 8], index % 8, 1);
}

void Baps::gen_w_baps() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint8_t> dis(0, 1);

    w_.resize(lambda_bytes_);
    for (auto& wi : w_) {
        wi = dis(gen);
    }
}

void Baps::get_z_baps() {
    bool flag = 0;
    z_.resize(32);
    for (int i = 0; i < 2 * lambda_bytes_; i++) {
        if (GET_BIT(y_[i / 8], i % 8) || flag) {
            setbit(&z_[i / 8], i % 8, 1);
        }
    }
}

void Baps::get_x_z_baps() {
    x_z.resize(64);
    std::copy(z_.begin(), z_.end(), x_z.data());
    std::copy(attribute_.begin(), attribute_.end(),
              x_z.data() + 2 * lambda_bytes_);
    std::copy(w_.begin(), w_.end(), x_z.data() + 3 * lambda_bytes_);
}

void Baps::get_b_baps() {
    b_.resize(k_);
    std::copy(policy_.begin(), policy_.end(), b_.data());
    std::copy(attribute_.begin(), attribute_.end(),
              b_.data() + 2 * lambda_bytes_);
}

void Baps::gen_G_baps() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint8_t> dis(0, 1);
    for (auto& g_i : G1_) {
        for (int i = 0; i < k_; i++) {
            g_i[i] = dis(gen);
        }
    }
    for (auto& g_i : G2_) {
        for (auto& value : g_i) {
            value = dis(gen);
        }
    }
}

void Baps::calculate_t_baps() {
    int tmp;
    for (int index = 0; index < k_; index++) {
        tmp = 0;
        for (int i = 0; i < k_; i++) {
            tmp += G2_[index][i] * GET_BIT(b_[i / 8], i % 8);
            tmp +=
                G1_[index][i] * GET_BIT(b_[i / 8], i % 8) * GET_BIT(b_[0], 0);
        }
        t_[index] = tmp % 2;
    }
}

void Baps::get_constrain_prover_disclose(const field::GF2_128* v,
                                         const uint8_t* witness,
                                         field::GF2_128* A_0,
                                         field::GF2_128* A_1) {
    for (int index = 0; index < k_; index++) {
        for (int i = 0; i < k_; i++) {
            if (G2_[index][i] == 1) A_1[index] += v[i];
            if (G1_[index][i] == 1) {
                A_0[index] += v[0] * v[i];
                if (GET_BIT(b_[i / 8], i % 8)) A_1[index] += v[0];
                if (GET_BIT(b_[0], 0)) A_1[index] += v[i];
            }
        }
    }
}

void Baps::get_constrain_verifier_disclose(const field::GF2_128* q,
                                           field::GF2_128 delta,
                                           field::GF2_128* B) {
    for (int index = 0; index < k_; index++) {
        if (t_[index] == 1) B[index] -= delta * delta;
        for (int i = 0; i < k_; i++) {
            if (G2_[index][i] == 1) B[index] += q[i] * delta;
            if (G1_[index][i] == 1) B[index] += q[0] * q[i];
        }
    }
}

void Baps::hash_mu_baps(std::vector<uint8_t>& mu,
                        const std::vector<uint8_t>& msg) {
    H1_context_t h1_ctx;
    H1_init(&h1_ctx, lambda_);
    H1_update(&h1_ctx, msg.data(), msg.size());
    H1_update(&h1_ctx, attribute.tree_root.data(), attribute.tree_root.size());
    H1_update(&h1_ctx, policy.tree_root.data(), policy.tree_root.size());
    H1_final(&h1_ctx, mu.data(), 2 * lambda_bytes_);
}

void Baps::hash_mu_disclose_baps(std::vector<uint8_t>& mu) {
    H1_context_t h1_ctx;
    H1_init(&h1_ctx, lambda_);
    H1_update(&h1_ctx, attribute.tree_root.data(), attribute.tree_root.size());
    H1_update(&h1_ctx, policy.tree_root.data(), policy.tree_root.size());
    H1_final(&h1_ctx, mu.data(), 2 * lambda_bytes_);
}
void Baps::gen_rootkey_iv_baps(const std::vector<uint8_t>& mu,
                               std::vector<uint8_t>& rootkey,
                               std::vector<uint8_t>& iv) {
    H3_context_t h3_ctx;
    H3_init(&h3_ctx, lambda_);
    H3_update(&h3_ctx, attribute_.data(), lambda_bytes_);
    H3_update(&h3_ctx, policy_.data(), lambda_bytes_);
    H3_update(&h3_ctx, mu.data(), lambda_bytes_ * 2);
    H3_final(&h3_ctx, rootkey.data(), lambda_bytes_, iv.data());
}

void Baps::baps_commit(const std::vector<uint8_t>& msg,
                       const std::vector<uint8_t>& r_bytes, uint8_t* witness) {
    field::GF2_128 m_1, m_2, r, rain_output_1, tmp, rain_output_2;
    r.from_bytes(r_bytes.data());
    r.to_bytes(witness);
    witness += 16UL;
    m_1.from_bytes(msg.data());
    m_1.to_bytes(witness);
    witness += 16UL;
    if (msg.size() == 2 * lambda_bytes_) {
        m_2.from_bytes(msg.data() + lambda_bytes_);
    }
    rain(m_1, r, rain_output_1, witness, true);
    witness += 16UL * (NUM_SBOXES - 1);
    tmp = rain_output_1 + r;
    tmp.to_bytes(witness);
    witness += 16UL;
    if (msg.size() == 2 * lambda_bytes_) {
        m_2.to_bytes(witness);
        witness += 16UL;
        rain(m_2, tmp, rain_output_2, witness, true);
        witness += 16UL * (NUM_SBOXES - 1);
        tmp = rain_output_2 + tmp;
        tmp.to_bytes(witness);
    }
}

void Baps::baps_sign(int index_attribute, int index_policy,
                     const std::vector<uint8_t>& msg, baps_signature_t* sig,double& during_time) {
    auto start_time = std::chrono::high_resolution_clock::now();
    // attribute.sign(index_attribute, msg, sig->sig_attribute);
    // policy.sign(index_policy, msg, sig->sig_policy);

    const unsigned int ell =
        2352 + 2 * lambda_ * 3 + lambda_ * (6 + 11) +
        (3 + NUM_SBOXES) * log2(num_members_) * lambda_ + lambda_ +
        (3 + NUM_SBOXES) * log2(num_members_ * 2) * lambda_ + lambda_;
    const unsigned int muti_times = 757 + 2 * lambda_ + 12 +
                                    (1 + NUM_SBOXES) * log2(num_members_) +
                                    (1 + NUM_SBOXES) * log2(num_members_ * 2);
    const unsigned int ell_bytes = ell / 8;
    const unsigned int ell_hat = ell + lambda_ * 2 + UNIVERSAL_HASH_B_BITS;
    const unsigned int ell_hat_bytes = ell_hat / 8;

    std::vector<uint8_t> mu(2 * lambda_bytes_);
    hash_mu_baps(mu, msg);
    std::vector<uint8_t> rootkey(lambda_bytes_);
    sig->iv.resize(iv_size_);
    gen_rootkey_iv_baps(mu, rootkey, sig->iv);

    std::vector<uint8_t> hcom(lambda_bytes_ * 2);
    std::vector<vec_com_t> vecCom(params_.tau);
    std::vector<uint8_t> u(ell_hat_bytes);
    std::vector<uint8_t*> V(lambda_);
    V[0] = new uint8_t[lambda_ * ell_hat_bytes];
    for (unsigned int i = 1; i < lambda_; ++i) {
        V[i] = V[0] + i * ell_hat_bytes;
    }

    sig->c.resize((params_.tau - 1) * ell_hat_bytes);

    vole_commit(rootkey.data(), sig->iv.data(), ell_hat, &(params_),
                hcom.data(), vecCom.data(), sig->c.data(), u.data(), V.data());

    std::vector<uint8_t> chall_1(5 * lambda_bytes_ + 8);
    hash_challenge_1(mu, hcom, sig->c, sig->iv, chall_1, lambda_, ell,
                     params_.tau);

    sig->u_tilde.resize(lambda_bytes_ + UNIVERSAL_HASH_B);

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
    gen_witness_255(y_.data(), witness.data(), 7);

    gen_witness_255(x_z.data(), witness.data() + 800 / 8, 8);

    gen_witness_dot_product(witness.data() + 2352 / 8, policy_.data(),
                            y_.data());

    std::vector<uint8_t> r(lambda_bytes_);

    baps_commit(attribute_, r, witness.data() + (2352 + 2 * lambda_ * 3) / 8);
    baps_commit(policy_, r,
                witness.data() + (2352 + 2 * lambda_ * 3 + lambda_ * 6) / 8);

    attribute.gen_witness(
        witness.data() + (2352 + 2 * lambda_ * 3 + lambda_ * (6 + 11)) / 8,
        index_attribute);

    policy.gen_witness(
        witness.data() +
            (2352 + 2 * lambda_ * 3 + lambda_ * (6 + 11) +
             (3 + NUM_SBOXES) * int(log2(num_members_)) * 128 + 128) /
                8,
        index_policy);

    sig->d.resize(ell_bytes);

    xor_u8_array(witness.data(), u.data(), sig->d.data(), ell_bytes);

    std::vector<uint8_t> chall_2(3 * lambda_bytes_ + 8);
    hash_challenge_2(chall_2, chall_1, sig->u_tilde, h_v, sig->d, lambda_, ell);

    std::vector<field::GF2_128> v_gf_128_vec(ell_hat);
    std::vector<field::GF2_128> v_combined_gf_128_vec(
        (6 + 11) + (3 + NUM_SBOXES) * log2(num_members_) + 1 +
        (3 + NUM_SBOXES) * log2(num_members_ * 2) + 1);

    convert_vec_to_field(V.data(), v_gf_128_vec.data(), ell_hat, lambda_);
    gen_combined_field_vec(
        v_gf_128_vec.data() + 2352 + 2 * lambda_ * 3,
        v_combined_gf_128_vec.data(),
        lambda_ * (6 + 11) + (3 + NUM_SBOXES) * log2(num_members_) * 128 + 128 +
            (3 + NUM_SBOXES) * log2(num_members_ * 2) * 128 + 128);

    std::vector<field::GF2_128> A_0(muti_times);
    std::vector<field::GF2_128> A_1(muti_times);

    get_constrain_prover(v_gf_128_vec.data(), witness.data(), A_0.data(),
                         A_1.data(), 255, 7);

    get_constrain_prover(v_gf_128_vec.data() + 800, witness.data() + 800 / 8,
                         A_0.data() + 255, A_1.data() + 255, 502, 8);

    get_constrain_prover_dot_product(
        v_gf_128_vec.data() + 2336 + 16, witness.data() + (2336 + 16) / 8,
        A_0.data() + 757, A_1.data() + 757, 2 * lambda_);

    baps_commit_prove(witness.data() + (2352 + 2 * lambda_ * 3) / 8,
                      v_combined_gf_128_vec.data(),
                      v_gf_128_vec.data() + 2352 + 2 * lambda_ * 3,
                      A_0.data() + 757 + 2 * lambda_,
                      A_1.data() + 757 + 2 * lambda_);

    path_prove(
        witness.data() + (2352 + 2 * lambda_ * 3 + lambda_ * (6 + 11)) / 8,
        v_combined_gf_128_vec.data() + 17,
        v_gf_128_vec.data() + 2352 + 2 * lambda_ * 3 + lambda_ * (6 + 11),
        attribute.rain_msg_, A_0.data() + 757 + 2 * lambda_ + 12,
        A_1.data() + 757 + 2 * lambda_ + 12, log2(num_members_));

    path_prove(
        witness.data() +
            (2352 + 2 * lambda_ * 3 + lambda_ * (6 + 11) +
             (3 + NUM_SBOXES) * int(log2(num_members_)) * 128 + 128) /
                8,
        v_combined_gf_128_vec.data() +
            int(17 + (3 + NUM_SBOXES) * log2(num_members_) + 1),
        v_gf_128_vec.data() + 2352 + 2 * lambda_ * 3 + lambda_ * (6 + 11) +
            (3 + NUM_SBOXES) * int(log2(num_members_)) * lambda_ + lambda_,
        attribute.rain_msg_,
        A_0.data() +
            int(757 + 2 * lambda_ + 12 + (1 + NUM_SBOXES) * log2(num_members_)),
        A_1.data() +
            int(757 + 2 * lambda_ + 12 + (1 + NUM_SBOXES) * log2(num_members_)),
        int(log2(2 * num_members_)));

    field::GF2_128 u_star;
    field::GF2_128 v_star;
    u_star.from_bytes(u.data() + ell_bytes);
    v_star = combine_bf128_vec(v_gf_128_vec.data() + ell);
    field::GF2_128 A_0_tilde = zk_hash(chall_2, A_0, v_star);
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

    auto total_time = end_time - start_time;
    during_time =  std::chrono::duration<double, std::milli>(total_time).count();

    // std::cout << "sign total time is : "
    //           << std::chrono::duration<double, std::milli>(total_time).count()
    //           << " ms" << std::endl;

    delete[] V[0];
}

bool Baps::baps_verify(const std::vector<uint8_t>& msg, baps_signature_t* sig, double& during_time) {
    auto start_time = std::chrono::high_resolution_clock::now();

    // if (!attribute.verify(msg, sig->sig_attribute)) return false;
    // if (!policy.verify(msg, sig->sig_policy)) return false;
    const unsigned int ell =
        2352 + 2 * lambda_ * 3 + lambda_ * (6 + 11) +
        (3 + NUM_SBOXES) * log2(num_members_) * lambda_ + lambda_ +
        (3 + NUM_SBOXES) * log2(num_members_ * 2) * lambda_ + lambda_;
    const unsigned int muti_times = 757 + 2 * lambda_ + 12 +
                                    (1 + NUM_SBOXES) * log2(num_members_) +
                                    (1 + NUM_SBOXES) * log2(num_members_ * 2);
    const unsigned int ell_bytes = ell / 8;
    const unsigned int ell_hat = ell + lambda_ * 2 + UNIVERSAL_HASH_B_BITS;
    const unsigned int ell_hat_bytes = ell_hat / 8;

    std::vector<uint8_t> mu(2 * lambda_bytes_);
    hash_mu_baps(mu, msg);

    std::vector<uint8_t*> Q(lambda_);
    std::vector<uint8_t> hcom(lambda_bytes_ * 2);
    Q[0] = new uint8_t[lambda_ * ell_hat_bytes];
    for (unsigned int i = 1; i < lambda_; ++i) {
        Q[i] = Q[0] + i * ell_hat_bytes;
    }
    vole_reconstruct(sig->iv.data(), sig->chall_3.data(), sig->pdec.data(),
                     sig->com.data(), hcom.data(), Q.data(), ell_hat, &params_);

    std::vector<uint8_t> chall_1(5 * lambda_bytes_ + 8);
    hash_challenge_1(mu, hcom, sig->c, sig->iv, chall_1, lambda_, ell,
                     params_.tau);

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
    delete[] Dtilde[0];

    std::vector<uint8_t> chall_2(3 * lambda_bytes_ + 8);
    hash_challenge_2(chall_2, chall_1, sig->u_tilde, h_v, sig->d, lambda_, ell);

    // if (compare_uint8_arrays(chall_2, sig->chall_2))
    //     std::cout << "yes" << std::endl;
    // else
    //     std::cout << "no" << std::endl;

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
    std::vector<field::GF2_128> q_combined_gf_128_vec(
        17 + (3 + NUM_SBOXES) * log2(num_members_) + 1 +
        (3 + NUM_SBOXES) * log2(num_members_ * 2) + 1);

    convert_vec_to_field(Q_.data(), q_gf_128_vec.data(), ell_hat, lambda_);

    gen_combined_field_vec(
        q_gf_128_vec.data() + 2352 + 2 * lambda_ * 3,
        q_combined_gf_128_vec.data(),
        lambda_ * (6 + 11) + (3 + NUM_SBOXES) * log2(num_members_) * 128 + 128 +
            (3 + NUM_SBOXES) * log2(num_members_ * 2) * 128 + 128);

    field::GF2_128 delta_field;
    delta_field.from_bytes(sig->chall_3.data());

    std::vector<field::GF2_128> B(muti_times);
    get_constrain_verifier(q_gf_128_vec.data(), delta_field, B.data(), 255, 7);
    get_constrain_verifier(q_gf_128_vec.data() + 800, delta_field,
                           B.data() + 255, 502, 8);
    get_constrain_verifier_dot_product(q_gf_128_vec.data() + 2336 + 16,
                                       delta_field, B.data() + 757,
                                       2 * lambda_);
    baps_commit_verify(q_combined_gf_128_vec.data(),
                       q_gf_128_vec.data() + 2352 + 2 * lambda_ * 3,
                       delta_field, B.data() + 757 + 2 * lambda_);

    path_verify(q_combined_gf_128_vec.data() + 17,
                q_gf_128_vec.data() + 2352 + 2 * lambda_ * 3 + 17 * lambda_,
                delta_field, attribute.rain_msg_,
                B.data() + 757 + 2 * lambda_ + 12, log2(num_members_));

    path_verify(q_combined_gf_128_vec.data() + 17 +
                    (3 + NUM_SBOXES) * int(log2(num_members_)) + 1,
                q_gf_128_vec.data() + 2352 + 2 * lambda_ * 3 + 17 * lambda_ +
                    (3 + NUM_SBOXES) * int(log2(num_members_)) * lambda_ +
                    lambda_,
                delta_field, attribute.rain_msg_,
                B.data() + 757 + 2 * lambda_ + 12 +
                    (1 + NUM_SBOXES) * int(log2(num_members_)),
                log2(2 * num_members_));

    // for (int i = 0; i < muti_times; i++) {
    //     if (B[i] != sig->A_0[i] + sig->A_1[i] * delta_field) {
    //         // std::cout << i << "fail" << std::endl;
    //     } else {
    //         std::cout << i << "pass" << std::endl;
    //     }
    // }
    field::GF2_128 zero;
    field::GF2_128 q_star;
    q_star = combine_bf128_vec(q_gf_128_vec.data() + ell);
    field::GF2_128 B_tilde = zk_hash(chall_2, B, q_star);

    field::GF2_128 A_0_tilde, A_1_tilde;
    std::vector<uint8_t> A_0_tilde_bytes(lambda_bytes_);

    A_1_tilde.from_bytes(sig->A_1_tilde_bytes.data());
    A_0_tilde = B_tilde - delta_field * A_1_tilde;
    A_0_tilde.to_bytes(A_0_tilde_bytes.data());

    std::vector<uint8_t> chall_3(lambda_bytes_);

    hash_challenge_3(chall_3, chall_2, sig->A_1_tilde_bytes, A_0_tilde_bytes,
                     lambda_);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto total_time = end_time - start_time;
    during_time =  std::chrono::duration<double, std::milli>(total_time).count();

    // std::cout << "verify total time is : "
    //           << std::chrono::duration<double, std::milli>(total_time).count()
    //           << " ms" << std::endl;

    delete[] Q[0];
    delete[] Q_[0];

    return memcmp(chall_3.data(), sig->chall_3.data(), lambda_bytes_) == 0;
}

void Baps::baps_disclose(baps_attestation_t* sig,double& during_time) {
    const unsigned int ell = k_;
    const unsigned int muti_times = k_;
    const unsigned int ell_bytes = ell / 8;
    const unsigned int ell_hat = ell + lambda_ * 2 + UNIVERSAL_HASH_B_BITS;
    const unsigned int ell_hat_bytes = ell_hat / 8;

    auto start_time = std::chrono::high_resolution_clock::now();
    std::vector<uint8_t> mu(2 * lambda_bytes_);
    hash_mu_disclose_baps(mu);
    std::vector<uint8_t> rootkey(lambda_bytes_);

    sig->iv.resize(iv_size_);
    gen_rootkey_iv_baps(mu, rootkey, sig->iv);

    std::vector<uint8_t> hcom(lambda_bytes_ * 2);
    std::vector<vec_com_t> vecCom(params_.tau);
    std::vector<uint8_t> u(ell_hat_bytes);
    std::vector<uint8_t*> V(lambda_);
    V[0] = new uint8_t[lambda_ * ell_hat_bytes];
    for (unsigned int i = 1; i < lambda_; ++i) {
        V[i] = V[0] + i * ell_hat_bytes;
    }

    sig->c.resize((params_.tau - 1) * ell_hat_bytes);

    vole_commit(rootkey.data(), sig->iv.data(), ell_hat, &(params_),
                hcom.data(), vecCom.data(), sig->c.data(), u.data(), V.data());

    std::vector<uint8_t> chall_1(5 * lambda_bytes_ + 8);
    hash_challenge_1(mu, hcom, sig->c, sig->iv, chall_1, lambda_, ell,
                     params_.tau);

    sig->u_tilde.resize(lambda_bytes_ + UNIVERSAL_HASH_B);

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
    witness = b_;

    sig->d.resize(ell_bytes);

    xor_u8_array(witness.data(), u.data(), sig->d.data(), ell_bytes);

    std::vector<uint8_t> chall_2(3 * lambda_bytes_ + 8);
    hash_challenge_2(chall_2, chall_1, sig->u_tilde, h_v, sig->d, lambda_, ell);

    std::vector<field::GF2_128> v_gf_128_vec(ell_hat);

    convert_vec_to_field(V.data(), v_gf_128_vec.data(), ell_hat, lambda_);

    std::vector<field::GF2_128> A_0(muti_times);
    std::vector<field::GF2_128> A_1(muti_times);
    get_constrain_prover_disclose(v_gf_128_vec.data(), witness.data(),
                                  A_0.data(), A_1.data());

    field::GF2_128 u_star;
    field::GF2_128 v_star;
    u_star.from_bytes(u.data() + ell_bytes);
    v_star = combine_bf128_vec(v_gf_128_vec.data() + ell);
    field::GF2_128 A_0_tilde = zk_hash(chall_2, A_0, v_star);
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

    auto total_time = end_time - start_time;
    during_time =  std::chrono::duration<double, std::milli>(total_time).count();

    // std::cout << "disclose total time is : "
    //           << std::chrono::duration<double, std::milli>(total_time).count()
    //           << " ms" << std::endl;

    delete[] V[0];
}

bool Baps::baps_disclose_verify(baps_attestation_t* sig,double& during_time) {
    const unsigned int ell = k_;
    const unsigned int muti_times = k_;
    const unsigned int ell_bytes = ell / 8;
    const unsigned int ell_hat = ell + lambda_ * 2 + UNIVERSAL_HASH_B_BITS;
    const unsigned int ell_hat_bytes = ell_hat / 8;

    auto start_time = std::chrono::high_resolution_clock::now();

    std::vector<uint8_t> mu(2 * lambda_bytes_);
    hash_mu_disclose_baps(mu);

    std::vector<uint8_t*> Q(lambda_);
    std::vector<uint8_t> hcom(lambda_bytes_ * 2);
    Q[0] = new uint8_t[lambda_ * ell_hat_bytes];
    for (unsigned int i = 1; i < lambda_; ++i) {
        Q[i] = Q[0] + i * ell_hat_bytes;
    }
    vole_reconstruct(sig->iv.data(), sig->chall_3.data(), sig->pdec.data(),
                     sig->com.data(), hcom.data(), Q.data(), ell_hat, &params_);

    std::vector<uint8_t> chall_1(5 * lambda_bytes_ + 8);
    hash_challenge_1(mu, hcom, sig->c, sig->iv, chall_1, lambda_, ell,
                     params_.tau);

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
    delete[] Dtilde[0];

    std::vector<uint8_t> chall_2(3 * lambda_bytes_ + 8);
    hash_challenge_2(chall_2, chall_1, sig->u_tilde, h_v, sig->d, lambda_, ell);

    // if (compare_uint8_arrays(chall_2, sig->chall_2))
    //     std::cout << "yes" << std::endl;
    // else
    //     std::cout << "no" << std::endl;

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
    convert_vec_to_field(Q_.data(), q_gf_128_vec.data(), ell_hat, lambda_);

    field::GF2_128 delta_field;
    delta_field.from_bytes(sig->chall_3.data());

    std::vector<field::GF2_128> B(muti_times);
    get_constrain_verifier_disclose(q_gf_128_vec.data(), delta_field, B.data());

    // for (int i = 0; i < muti_times; i++) {
    //     if (B[i] != sig->A_0[i] + sig->A_1[i] * delta_field) {
    //         std::cout << i << "fail" << std::endl;
    //     } else {
    //         std::cout << i << "pass" << std::endl;
    //     }
    // }

    field::GF2_128 zero;
    field::GF2_128 q_star;
    q_star = combine_bf128_vec(q_gf_128_vec.data() + ell);
    field::GF2_128 B_tilde = zk_hash(chall_2, B, q_star);

    field::GF2_128 A_0_tilde, A_1_tilde;
    std::vector<uint8_t> A_0_tilde_bytes(lambda_bytes_);

    A_1_tilde.from_bytes(sig->A_1_tilde_bytes.data());
    A_0_tilde = B_tilde - delta_field * A_1_tilde;
    A_0_tilde.to_bytes(A_0_tilde_bytes.data());

    std::vector<uint8_t> chall_3(lambda_bytes_);

    hash_challenge_3(chall_3, chall_2, sig->A_1_tilde_bytes, A_0_tilde_bytes,
                     lambda_);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto total_time = end_time - start_time;
    during_time =  std::chrono::duration<double, std::milli>(total_time).count();

    // std::cout << "disclose verify total time is : "
    //           << std::chrono::duration<double, std::milli>(total_time).count()
    //           << " ms" << std::endl;

    delete[] Q[0];
    delete[] Q_[0];

    return memcmp(chall_3.data(), sig->chall_3.data(), lambda_bytes_) == 0;
}