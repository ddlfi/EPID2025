#include "signature.h"

#include <chrono>
#include <cmath>
#include <random>

void Signature::calculate_params() {
    int tmp = m_;
    while (tmp > 1) {
        ell_ += tmp + tmp / 2;
        muti_num_ += tmp / 2;
        tmp /= 2;
    }
    ell_++;
    ell_ += (8 - ell_ % 8);
}

void Signature::gen_matrix() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint8_t> dis(0, 255);
    for (auto& m_i : matrix_) {
        for (auto& m_i_j : m_i) m_i_j = dis(gen);
    }
}

void Signature::gen_x() {
    const int num_bytes = (m_ + 7) / 8;
    x_.resize(num_bytes);

    std::vector<int> positions(m_);
    std::iota(positions.begin(), positions.end(), 0);

    std::random_device rd;
    std::mt19937 rng(rd());
    std::shuffle(positions.begin(), positions.end(), rng);

    for (int i = 0; i < hamming_weight; i++) {
        const int pos = positions[i];
        setbit(&x_[pos / 8], pos % 8,1);
    }
}

void Signature::calculate_y() {

    for (const auto& row : matrix_) {
        uint8_t acc = 0;
        for (size_t byte_idx = 0; byte_idx < (m_ + 7) / 8; byte_idx++) {
            acc ^= (row[byte_idx] & x_[byte_idx]);
        }
        y_.push_back(parity(acc));
    }
}

void Signature::hash_mu(const std::vector<uint8_t>& msg,
                        std::vector<uint8_t>& mu) {
    H1_context_t h1_ctx;
    H1_init(&h1_ctx, lambda_);
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
                               std::vector<uint8_t>& rootkey,
                               std::vector<uint8_t>& iv) {
    H3_context_t h3_ctx;
    H3_init(&h3_ctx, lambda_);
    H3_update(&h3_ctx, mu.data(), lambda_bytes_ * 2);
    H3_final(&h3_ctx, rootkey.data(), lambda_bytes_, iv.data());
}

void Signature::sign(const std::vector<uint8_t>& msg, signature_t* sig) {
    const unsigned int ell = ell_;
    const unsigned int muti_times = muti_num_;
    const unsigned int ell_bytes = ell / 8;
    const unsigned int ell_hat = ell + lambda_ * 2 + UNIVERSAL_HASH_B_BITS;
    const unsigned int ell_hat_bytes = ell_hat / 8;

    auto start_time = std::chrono::high_resolution_clock::now();
    std::vector<uint8_t> mu(2 * lambda_bytes_);
    hash_mu(msg, mu);

    std::vector<uint8_t> rootkey(lambda_bytes_);
    sig->iv.resize(iv_size_);
    gen_rootkey_iv(mu, rootkey, sig->iv);

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
    hash_challenge_1(mu, hcom, sig->c, sig->iv, chall_1, ell, params_.tau);

    sig->u_tilde.resize(lambda_bytes_ + UNIVERSAL_HASH_B);
    vole_hash(sig->u_tilde.data(), chall_1.data(), u.data(), ell, lambda_);

    std::vector<uint8_t> h_v(lambda_bytes_ * 2);
    {
        H1_context_t h1_ctx_1;
        H1_init(&h1_ctx_1, lambda_);
        std::vector<uint8_t> V_tilde(lambda_bytes_ + UNIVERSAL_HASH_B);
        for (unsigned int i = 0; i != lambda_; ++i) {
            vole_hash(V_tilde.data(), chall_1.data(), V[i], ell, lambda_);
            H1_update(&h1_ctx_1, V_tilde.data(),
                      lambda_bytes_ + UNIVERSAL_HASH_B);
        }
        H1_final(&h1_ctx_1, h_v.data(), lambda_bytes_ * 2);
    }

    sig->d.resize(ell_bytes);
    std::vector<uint8_t> witness(ell_bytes);
    gen_witness(x_.data(), witness.data(), m_);
    xor_u8_array(witness.data(), u.data(), sig->d.data(), ell_bytes);

    std::vector<uint8_t> chall_2(3 * lambda_bytes_ + 8);
    hash_challenge_2(chall_2, chall_1, sig->u_tilde, h_v, sig->d, lambda_, ell);

    std::vector<field::GF2_128> v_gf_128_vec(ell_hat);
    convert_vec_to_field(V.data(), v_gf_128_vec.data(), ell_hat, lambda_);

    sig->A_linear.resize(k_);
    get_constrain_prover_linear(matrix_,v_gf_128_vec.data(),sig->A_linear.data(),m_);

    std::vector<field::GF2_128> A_0(muti_times);
    std::vector<field::GF2_128> A_1(muti_times);
    sig->hamming_weight_vec.resize(std::floor(log2(m_) + 1e-10));
    sig->hamming_weight_v.resize(std::floor(log2(m_) + 1e-10));

    get_hamming_weight(witness.data(), v_gf_128_vec.data(),
                       sig->hamming_weight_vec.data(),
                       sig->hamming_weight_v.data(), m_);

    get_constrain_prover(v_gf_128_vec.data(), witness.data(), A_0.data(),
                         A_1.data(), muti_times, m_);

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

    sig->pdec.resize(params_.tau);
    sig->com.resize(params_.tau);

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
    auto during_time =
        std::chrono::duration<double, std::milli>(total_time).count();
    std::cout << "sign total time is : "
              << std::chrono::duration<double, std::milli>(total_time).count()
              << " ms" << std::endl;

    delete[] V[0];
}

bool Signature::verify(const std::vector<uint8_t>& msg,
                       const signature_t* sig) {
    const unsigned int ell = ell_;
    const unsigned int muti_times = muti_num_;
    const unsigned int ell_bytes = ell / 8;
    const unsigned int ell_hat = ell + lambda_ * 2 + UNIVERSAL_HASH_B_BITS;
    const unsigned int ell_hat_bytes = ell_hat / 8;

    std::vector<uint8_t> mu(2 * lambda_bytes_);
    hash_mu(msg, mu);

    std::vector<uint8_t*> Q(lambda_);
    std::vector<uint8_t> hcom(lambda_bytes_ * 2);
    Q[0] = new uint8_t[lambda_ * ell_hat_bytes];
    for (unsigned int i = 1; i < lambda_; ++i) {
        Q[i] = Q[0] + i * ell_hat_bytes;
    }
    vole_reconstruct(sig->iv.data(), sig->chall_3.data(), sig->pdec.data(),
                     sig->com.data(), hcom.data(), Q.data(), ell_hat, &params_);

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
    delete[] Dtilde[0];

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
    convert_vec_to_field(Q_.data(), q_gf_128_vec.data(), ell_hat, lambda_);

    field::GF2_128 delta_field;
    delta_field.from_bytes(sig->chall_3.data());
    

    std::vector<field::GF2_128> B_linear(k_);
    get_constrain_verifier_linear(matrix_,y_,q_gf_128_vec.data(),delta_field,B_linear.data(),m_);
    if(!verify_constrain_linear(sig->A_linear.data(),B_linear.data(),k_)){
        return false;
    }

    // check hamming weight
    if (!verify_hamming_weight(hamming_weight, sig->hamming_weight_vec.data(),
                               sig->hamming_weight_v.data(),
                               q_gf_128_vec.data(), delta_field, m_)) {
        return false;
    }



    std::vector<field::GF2_128> B(muti_times);
    get_constrain_verifier(q_gf_128_vec.data(), delta_field, B.data(),
                           muti_times, m_);

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

    delete[] Q[0];
    delete[] Q_[0];

    return memcmp(chall_3.data(), sig->chall_3.data(), lambda_bytes_) == 0;
}
