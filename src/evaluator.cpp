#include "fHE/evaluator.hpp"
#include <array>
#include "fHE/base_converter.hpp"
#include "fHE/basis_extension.hpp"
#include "fHE/cipher.hpp"
#include "fHE/encoder.hpp"
#include "fHE/mult_key.hpp"
#include "fHE/prime_bundles.hpp"
#include "fHE/rot_key.hpp"
#include "fHE/special_prime_chain.hpp"
#include "fHE/yell.hpp"

namespace fHE {

struct Evaluator::Impl {
    // BaseConverter bconv;

    std::shared_ptr<PrimeBundles> bundles_;
    std::shared_ptr<SpecialPrimeChain> chain_;
    BasisExtension bext_;

    explicit Impl(std::shared_ptr<PrimeBundles> bundles, std::shared_ptr<SpecialPrimeChain> chain)
        : bundles_(bundles), chain_(chain), bext_(bundles, chain) {}

    inline bool can_add_with(Cipher const &c0, Cipher const &c1) const {
        // TODO check scale
        if (c0.moduli_count() != c1.moduli_count()) return false;
        return true;
    }

    inline bool can_mul_with(Cipher const &c0, Cipher const &c1) const {
        if (c0.moduli_count() != c1.moduli_count()) return false;
        if (!(c0.canonical() && c1.canonical())) return false;
        return true;
    }

    bool add_inplace(Cipher *rop, Cipher const &c1) const;

    bool multiply_inplace(Cipher *rop, Cipher const &c1) const;

    bool relinerize_inplace(Cipher *rop, MultKey const &mkey) const;

    bool rotate_slots_inplace(Cipher *rop, RotKey const &rkey) const;

    bool rotate_slots_inplace_non_ntt(Cipher *rop, RotKey const &rkey) const;

private:
    struct Mul {
        void cross_mult(context::poly_t ans[3], Cipher const &c0, Cipher const &c1) const;

        void reLinearize(std::array<context::poly_t *, 2> &rop,
                         std::array<const context::poly_t *, 3> ct,
                         MultKey const &evk) const;

        // void rescale_ntt(std::array<context::poly_t *, 2> rop,
        //                  std::array<context::poly_t *, 2> ct,
        //                  BaseConverter const& bconv) const;
    };

    void apply_rotation_key(std::array<context::poly_t *, 2> rop, context::poly_t const &d, RotKey const &rotk) const;

    void apply_rotation_key_non_ntt(std::array<context::poly_t *, 2> rop,
                                    std::array<context::poly_t *, 2> d,
                                    RotKey const &rotk) const;
};

Evaluator::Evaluator(std::shared_ptr<Encoder> encoder,
                     std::shared_ptr<PrimeBundles> bundles,
                     std::shared_ptr<SpecialPrimeChain> chain)
    : impl_(std::make_shared<Impl>(bundles, chain)), encoder(encoder) {}

Evaluator::~Evaluator() {}

bool Evaluator::add(Cipher *rop, Cipher const &c1) const {
    if (!impl_) return false;
    return impl_->add_inplace(rop, c1);
}

bool Evaluator::multiply(Cipher *rop, Cipher const &c1) const {
    if (!impl_) return false;
    return impl_->multiply_inplace(rop, c1);
}

bool Evaluator::multiply_plain(Cipher *rop, context::poly_t const &p, const double plain_scale) const {
    if (!rop) return false;
    size_t nm = p.moduli_count();
    if (rop->moduli_count() != nm) return false;
    yell::ops::mulmod mulmod;
    for (size_t cm = 0; cm < rop->moduli_count(); ++cm) {
        auto op0 = rop->bx.ptr_at(cm);
        auto op1 = rop->ax.ptr_at(cm);
        auto op2 = p.cptr_at(cm);
        for (size_t d = 0; d < context::degree; ++d) {
            auto _op2 = *op2++;
            mulmod.compute(*op0++, _op2, cm);
            mulmod.compute(*op1++, _op2, cm);
        }
    }

    if (plain_scale > 0.) rop->scale(rop->scale() * plain_scale);
    return true;
}

bool Evaluator::relinearize(Cipher *rop, MultKey const &mkey) const {
    if (!impl_) return false;
    return impl_->relinerize_inplace(rop, mkey);
}

bool Evaluator::multiply_and_relin(Cipher *rop, Cipher const &c1, MultKey const &mkey) const {
    if (!impl_) return false;
    if (impl_->multiply_inplace(rop, c1))
        return impl_->relinerize_inplace(rop, mkey);
    else
        return false;
}

bool Evaluator::rotate_slots(Cipher *rop, RotKey const &rkey) const {
    if (!impl_ || !rop) return false;
    if (!rop->canonical()) {
        std::cerr << "Can only rotate canonical ciphertext, call relinearize first." << std::endl;
        return false;
    }
    return impl_->rotate_slots_inplace(rop, rkey);
}

bool Evaluator::rotate_slots_non_ntt(Cipher *rop, RotKey const &rkey) const {
    if (!impl_ || !rop) return false;
    if (!rop->canonical()) {
        std::cerr << "Can only rotate canonical ciphertext, call relinearize first." << std::endl;
        return false;
    }
    return impl_->rotate_slots_inplace_non_ntt(rop, rkey);
}

bool Evaluator::rotate_slots(Cipher *rop, int offset, RotKeySet<RotKey> const &rkeys) const {
    if (!impl_ || !rop) return false;
    if (!rop->canonical()) {
        std::cerr << "Can only rotate canonical ciphertext, call relinearize first." << std::endl;
        return false;
    }
    constexpr size_t nr_slots = context::degree >> 1u;
    bool sign                 = offset < 0;
    offset                    = std::abs(offset) % nr_slots;
    if (sign) offset = nr_slots - offset;
    int steps = (int)(std::log2((double)offset) + 1);
    for (int i = 0; i < steps; ++i) {
        if (offset & (1 << i)) {
            auto rotkey = rkeys.get(1 << i);
            if (!rotkey) {
                std::cerr << "No such rotation key for offset " << (1 << i) << "\n";
                return false;
            }
            rotate_slots(rop, *rotkey);
        }
    }
    return true;
}

template <class RotKey>
bool Evaluator::rotate_slots_non_ntt(Cipher *rop, int offset, RotKeySet<RotKey> const &rkeys) const {
    if (!impl_ || !rop) return false;
    if (!rop->canonical()) {
        std::cerr << "Can only rotate canonical ciphertext, call relinearize first." << std::endl;
        return false;
    }
    constexpr size_t nr_slots = context::degree >> 1u;
    offset                    = std::abs(offset) % nr_slots;
    const int steps           = (int)(std::log2((double)offset) + 1);
    for (int i = 0; i < steps; ++i) {
        if (offset & (1 << i)) {
            RotKey const *rotkey = rkeys.get(1 << i);
            if (!rotkey) {
                std::cerr << "No such rotation key for offset " << (1 << i) << "\n";
                return false;
            }
            rotate_slots_non_ntt(rop, *rotkey);
        }
    }
    return true;
}

template <class RotKey>
bool Evaluator::replicate(Cipher *rop, size_t pos, RotKeySet<RotKey> const &rkeys) const {
    if (!impl_ || !rop || !encoder) return false;
    if (!rop->canonical()) {
        std::cerr << "Can only replicate canonical ciphertext, call relinearize first." << std::endl;
        return false;
    }
    constexpr size_t nr_slots = context::degree >> 1;
    constexpr size_t logn     = yell::static_log2<nr_slots>::value;
    //! replicate uses one masking and logN rotations and additions.
    pos %= (nr_slots);
    std::vector<double> mask(nr_slots, 0.);
    mask[pos] = 1.;

    double scale = context::encoder_scale / 8.;
    context::poly_t mask_poly(rop->moduli_count());
    encoder->encode(&mask_poly, mask, scale);
    mask_poly.forward_lazy();
    multiply_plain(rop, mask_poly, scale);
    //! rotate above the power-basis is faster
    rop->bx.backward();
    rop->ax.backward();

    for (size_t i = 0; i < logn; ++i) {
        fHE::Cipher tmp(*rop);
        rotate_slots_non_ntt(&tmp, 1 << i, rkeys);
        add(rop, tmp);  // add is also ok on the power-basis.
    }

    rop->bx.forward();
    rop->ax.forward();
    return true;
}

// Implementation details go from here
bool Evaluator::Impl::add_inplace(Cipher *rop, Cipher const &c1) const {
    if (!rop) return false;
    if (!can_add_with(*rop, c1)) {
        std::cerr << "Can not add_inplace." << std::endl;
        return false;
    }

    if (!rop->canonical() && c1.canonical()) {
        std::cerr << "Can not add_inplace" << std::endl;
        return false;
    }
    rop->bx += c1.bx;
    rop->ax += c1.ax;

    if (!c1.canonical()) {
        if (rop->canonical())
            rop->mul_aux = std::make_shared<context::poly_t>(*c1.mul_aux);
        else
            (*rop->mul_aux) += (*c1.mul_aux);
    }

    if (rop->scale() == 1.) rop->scale(c1.scale());
    return true;
}

bool Evaluator::Impl::rotate_slots_inplace(Cipher *rop, RotKey const &rkey) const {
    if (!rop) return false;
    if (rkey.galois() == 0) return true;

    apply_galois(&(rop->bx), rkey.galois());
    apply_galois(&(rop->ax), rkey.galois());

    context::poly_t aux(rop->ax);
    std::array<context::poly_t *, 2> result {&rop->bx, &rop->ax};
    result[1]->clear();

    apply_rotation_key(result, aux, rkey);
    return true;
}

bool Evaluator::Impl::rotate_slots_inplace_non_ntt(Cipher *rop, RotKey const &rkey) const {
    if (!rop) return false;
    if (rkey.galois() == 0) return true;
    apply_galois_non_ntt(&(rop->bx), rkey.galois());
    apply_galois_non_ntt(&(rop->ax), rkey.galois());
    std::array<context::poly_t *, 2> result {&rop->bx, &rop->ax};
    std::array<context::poly_t *, 2> input {&rop->bx, &rop->ax};
    apply_rotation_key_non_ntt(result, input, rkey);
    return true;
}

//! Input polynomials [d] should be in the ntt domain.
//! Resulting polynomials [rop] are also in the ntt domain.

static inline size_t div_ceil(size_t a, size_t b) {
    assert(b > 0);
    return (a + b - 1) / b;
}

/// (out[0] + rotk.beta * aux, out[1] + rotk.alpha * aux)
void Evaluator::Impl::apply_rotation_key(std::array<context::poly_t *, 2> out,
                                         const context::poly_t &aux,
                                         RotKey const &rotk) const {
    if (!out[0] || !out[1]) throw std::invalid_argument("Nullptr error.");
    using poly_t = context::poly_t;
    using T      = yell::params::value_type;
    using gT     = yell::params::gt_value_type;

    const size_t max_n_nrml_primes = bundles_->max_n_nrml_primes();
    const size_t cur_n_nrml_primes = aux.moduli_count();  // the current number of normal primes
    const size_t degree            = aux.degree;
    const size_t max_n_bundles     = bundles_->n_bundles();
    const size_t bundle_size       = bundles_->bundle_size();
    const size_t n_spcl_primes     = bundle_size;
    const size_t cur_n_bundles     = div_ceil(cur_n_nrml_primes, bundle_size);  // the curent number of prime bundles

    assert(out[0]->moduli_count() == cur_n_nrml_primes);
    assert(out[1]->moduli_count() == cur_n_nrml_primes);

    /// RNS-decompose
    /// `rns_decomps` is `cur_n_bundles` polynomails, and each is using `bundle_size` moduli.
    /// Invariant: cur_n_bundles * bundle_size = cur_n_nrml_primes.
    std::vector<poly_t> rns_decomps;
    for (size_t j = 0; j < cur_n_bundles; ++j) {
        const size_t n_moduli = std::min(bundle_size * (j + 1), cur_n_nrml_primes) - bundle_size * j;
        rns_decomps.emplace_back(n_moduli);
        poly_t &rns_decomp = rns_decomps.back();  // point to the fresh poly.
        for (size_t i = 0; i < n_moduli; ++i) {
            size_t nrml_prime_idx = j * bundle_size + i;
            assert(nrml_prime_idx < cur_n_nrml_primes);
            yell::ops::mulmod mulmod;
            /// the j-th puncture bundle product under the 'nrml_prime_idx' prime.
            T punch_prod     = bundles_->puncture_product({.puncture_idx = j, .nrml_prime_idx = nrml_prime_idx});
            T punch_prod_inv = yell::math::inv_mod_prime(punch_prod, nrml_prime_idx);
            assert(punch_prod_inv > 0);

            if (punch_prod_inv > 1) {
                /// multiply 'remain_prod' to cmult[2] using shoup's trick for faster multiplication.
                T shoup = yell::ops::shoupify(punch_prod_inv, nrml_prime_idx);
                yell::ops::mulmod_shoup mulmod_s;
                std::transform(aux.cptr_at(nrml_prime_idx), aux.cptr_end(nrml_prime_idx), rns_decomp.ptr_at(i),
                               [&mulmod_s, &punch_prod_inv, &shoup, &nrml_prime_idx](T v) {
                                   return mulmod_s(v, punch_prod_inv, shoup, nrml_prime_idx);
                               });
            } else {
                /// Since punch_prod_inv == 1, we just copy cmult[2].
                size_t nbytes = sizeof(T) * degree;
                std::memcpy(rns_decomp.ptr_at(i), aux.cptr_at(nrml_prime_idx), nbytes);
            }
            /// back to power-basis for the following Mod-up operations.
            yell::ntt<context::degree>::backward(rns_decomp.ptr_at(i), nrml_prime_idx);
        }
    }

    /// Mod-up operations here. For each set of bundle primes, adding special primes.
    std::vector<poly_t> rns_decomps_mod_up;
    rns_decomps_mod_up.reserve(cur_n_bundles);
    const size_t n_total_moduli = n_spcl_primes + cur_n_nrml_primes;
    for (size_t j = 0; j < rns_decomps.size(); ++j) {
        /// the j-th set of bundle primes goes through
        /// p_{j * bundle_size + 0}, p_{j * bundle_size + 1}, ..., p_{(j+1) * bundle_size - 1}
        rns_decomps_mod_up.emplace_back(n_total_moduli);  // [noraml_part || special_part]
        auto dst_poly_ptr = &rns_decomps_mod_up.back();
        BasisExtension::ModupParams parms {
            .bundle_idx = j, .n_nrml_primes = cur_n_nrml_primes, .n_spcl_primes = n_spcl_primes};
        bext_.mod_up(dst_poly_ptr, rns_decomps[j], parms);

        /// Forward to NTT-basis for multiplication.
        for (size_t k = 0; k < n_total_moduli; ++k) {
            size_t prime_idx = k < cur_n_nrml_primes ? k : max_n_nrml_primes + k - cur_n_nrml_primes;
            yell::ntt<context::degree>::forward(dst_poly_ptr->ptr_at(k), prime_idx);
        }
    }

    /// Inner product.
    poly_t cmult_evk[2] {poly_t(n_total_moduli), poly_t(n_total_moduli)};
    std::vector<gT> lazy_mult[2] {std::vector<gT>(degree), std::vector<gT>(degree)};
    for (size_t k = 0; k < n_total_moduli; ++k) {
        bool is_spcl_prime = k >= cur_n_nrml_primes;
        size_t prime_idx   = is_spcl_prime ? max_n_nrml_primes + k - cur_n_nrml_primes : k;
        std::memset(lazy_mult[0].data(), 0, sizeof(gT) * degree);
        std::memset(lazy_mult[1].data(), 0, sizeof(gT) * degree);
        for (size_t j = 0; j < cur_n_bundles; ++j) {
            const T *ctx_ptr   = rns_decomps_mod_up[j].cptr_at(k);
            const T *beta_ptr  = rotk.beta_at(j).cptr_at(prime_idx);
            const T *alpha_ptr = rotk.alpha_at(j).cptr_at(prime_idx);
            gT *lazy_mult0_ptr = lazy_mult[0].data();
            gT *lazy_mult1_ptr = lazy_mult[1].data();

            for (size_t d = 0; d < degree; ++d, ++ctx_ptr) {
                *lazy_mult0_ptr++ += (gT)(*beta_ptr++) * (*ctx_ptr);
                *lazy_mult1_ptr++ += (gT)(*alpha_ptr++) * (*ctx_ptr);
            }
        }

        /// lazy reduction
        for (auto b : {0, 1}) {
            std::transform(lazy_mult[b].cbegin(), lazy_mult[b].cend(), cmult_evk[b].ptr_at(k), [prime_idx](gT v) -> T {
                yell::ops::barret_reduction(&v, prime_idx);
                return (T)v;
            });
            if (is_spcl_prime) {
                /// Backward to power-basis for mod-down. Note that, we only convert the special prime part.
                /// The following bext.mod_down() operation will take care of the normal part.
                yell::ntt<context::degree>::backward(cmult_evk[b].ptr_at(k), prime_idx);
            }
        }
    }

    /// Mod-down & add to ctxt.
    for (size_t b : {0, 1}) {
        bext_.mod_down_inplace(&cmult_evk[b], {.n_nrml_primes = cur_n_nrml_primes, .n_spcl_primes = n_spcl_primes});
        yell::ops::addmod addmod;
        /// Note: because the number of moduli of `out` is smaller than `cmult_evk`
        /// We only take cares of the normal prime part.
        for (size_t cm = 0; cm < cur_n_nrml_primes; ++cm) {
            std::transform(out[b]->cptr_at(cm), out[b]->cptr_end(cm), cmult_evk[b].cptr_at(cm), out[b]->ptr_at(cm),
                           [cm, addmod](const T a, const T b) -> T { return addmod(a, b, cm); });
        }
    }
}

//! Input polynomials [d] should be in the power domain.
//! Resulting polynomials [rop] are also in the power domain.
void Evaluator::Impl::apply_rotation_key_non_ntt(std::array<context::poly_t *, 2> rop,
                                                 std::array<context::poly_t *, 2> d,
                                                 RotKey const &rotk) const {
    // constexpr size_t L = context::nr_ctxt_moduli;
    // constexpr size_t K = context::nr_sp_primes;
    // if (!rop[0] || !rop[1] || !d[0] || !d[1])
    //   throw std::invalid_argument("Nullptr error.");
    // const size_t Li = rop[0]->moduli_count();
    // if (Li > L)
    //   throw std::invalid_argument("Too many moduli.");
    // assert(d[1]->moduli_count() == Li);
    // context::poly_t d_rotk[2] = { context::poly_t(Li + K),
    //                               context::poly_t(Li + K) };
    //
    // //! convert back to power-basis for the mod up operation.
    // bconv.approximated_mod_up(&d_rotk[0], (*d[1]));
    // //! forward to ntt domain for multiplication but now we have special moduli in d_rotk[0].
    // for (size_t cm = 0; cm < Li + K; ++cm) {
    //   const size_t moduli_idx = cm < Li ? cm : context::index_sp_prime(cm - Li);
    //   yell::ntt<context::degree>::forward(d_rotk[0].ptr_at(cm), moduli_idx);
    // }
    // d_rotk[1] = d_rotk[0]; // copy
    //
    // //! multiply the d[1] to rotation key over the extened moduli.
    // //! d[1] consists of Li (Li <= L) normal moduli and K special moduli,
    // //! however the rotation key consists of L normal moduli and K special moduli.
    // //! We need to carefully use the right moduli when doing the multiplication.
    // for (size_t cm = 0; cm < Li + K; ++cm) {
    //   const size_t moduli_idx = cm < Li ? cm : context::index_sp_prime(cm - Li);
    //   auto  beta = rotk.beta.cptr_at(moduli_idx);
    //   auto alpha = rotk.alpha.cptr_at(moduli_idx);
    //   auto op0   = d_rotk[0].ptr_at(cm);
    //   auto op1   = d_rotk[1].ptr_at(cm);
    //   yell::ops::mulmod mulmod;
    //   for (size_t d = 0; d < context::degree; ++d) {
    //     mulmod.compute(*op0++, *beta++, moduli_idx);
    //     mulmod.compute(*op1++, *alpha++, moduli_idx);
    //   }
    // }
    //
    // //! backward to power-basis for the mod down operation (with special moduli)
    // for (size_t cm = 0; cm < Li + K; ++cm) {
    //   const size_t moduli_idx = cm < Li ? cm : context::index_sp_prime(cm - Li);
    //   yell::ntt<context::degree>::backward(d_rotk[0].ptr_at(cm), moduli_idx);
    //   yell::ntt<context::degree>::backward(d_rotk[1].ptr_at(cm), moduli_idx);
    // }
    //
    // auto cpy_d0(*d[0]); // d[0], rop[0] might point to the same object,
    // //! mod down
    // for (int i : {0, 1})
    //   bconv.approximated_mod_down(rop[i], d_rotk[i]); // rop[i] is now in the power-basis
    // (*rop[0]) += cpy_d0;
}

bool Evaluator::Impl::multiply_inplace(Cipher *rop, Cipher const &c1) const {
    if (!rop) return false;
    if (!can_mul_with(*rop, c1)) {
        std::cerr << "Can not multiply_inplace" << std::endl;
        return false;
    }
    using poly_t           = context::poly_t;
    const size_t Li        = rop->moduli_count();
    const double new_scale = rop->scale() / yell::params::P[Li - 1] * c1.scale();
    poly_t ct[3] {poly_t(Li), poly_t(Li), poly_t(Li)};

    Mul mul;
    mul.cross_mult(ct, *rop, c1);
    rop->bx      = std::move(ct[0]);
    rop->ax      = std::move(ct[1]);
    rop->mul_aux = std::make_shared<context::poly_t>(ct[2]);
    rop->scale(new_scale);
    return true;
}

bool Evaluator::Impl::relinerize_inplace(Cipher *rop, MultKey const &mkey) const {
    if (!rop) return false;
    if (rop->canonical()) return false;
    using poly_t    = context::poly_t;
    const size_t Li = rop->moduli_count();
    //! Relinear
    //! NOTE: rop here are in the NTT domain.
    Mul mul;
    std::array<poly_t *, 2> tuple {&rop->bx, &rop->ax};
    std::array<const poly_t *, 3> triple {&rop->bx, &rop->ax, rop->mul_aux.get()};
    mul.reLinearize(tuple, triple, mkey);
    //! Results from reLinearize_via_bit_decompose are in the NTT domain.
    Cipher tmp(Li - 1);
    std::array<poly_t *, 2> res {&tmp.bx, &tmp.ax};
    // mul.rescale_ntt(res, tuple, bconv);
    rop->bx      = std::move(tmp.bx);
    rop->ax      = std::move(tmp.ax);
    rop->mul_aux = nullptr;
    return true;
}

void Evaluator::Impl::Mul::cross_mult(context::poly_t ans[3], Cipher const &c0, Cipher const &c1) const {
    using T              = yell::params::value_type;
    using gT             = yell::params::gt_value_type;
    const size_t nmoduli = ans[0].moduli_count();
    for (size_t cm = 0; cm < nmoduli; ++cm) {
        auto dst    = ans[1].ptr_at(cm);
        auto b0_ptr = c0.bx.cptr_at(cm);
        auto a0_ptr = c0.ax.cptr_at(cm);
        auto b1_ptr = c1.bx.cptr_at(cm);
        auto a1_ptr = c1.ax.cptr_at(cm);

        for (size_t d = 0; d < context::degree; ++d) {
            gT v0   = (gT)(*b0_ptr++) * (*a1_ptr++);
            auto v1 = (gT)(*b1_ptr++) * (*a0_ptr++);
            v0 += v1;
            yell::ops::barret_reduction(&v0, cm);
            *dst++ = (T)v0;
        }
    }

    ans[0] = c0.bx;
    ans[0] *= c1.bx;
    ans[2] = c0.ax;
    ans[2] *= c1.ax;
}

void Evaluator::Impl::Mul::reLinearize(std::array<context::poly_t *, 2> &rop,
                                       std::array<const context::poly_t *, 3> ct,
                                       MultKey const &evk) const {
    using T                 = yell::params::value_type;
    using gT                = yell::params::gt_value_type;
    constexpr size_t degree = context::degree;
    constexpr size_t bytes  = degree * sizeof(T);
    const size_t Li         = ct[2]->moduli_count();
    auto ct2_power(*ct[2]);
    ct2_power.backward();
    //! evk.beta * ct[2], and evk.alpha * ct[2]
    std::array<T, degree> ct2;
    for (size_t i = 0; i < Li; ++i) {
        //! compute the i-th moduli of evk * ct[2]
        std::array<gT, degree> evk_ct2[2] {};
        for (size_t j = 0; j < Li; ++j) {
            std::memcpy(ct2.data(), ct2_power.cptr_at(j), bytes);
            yell::ntt<degree>::forward(ct2.data(), i);
            //! Lazy reduction for evk.beta * ct[2]
            //! Lazy reduction for evk.alpha * ct[2]
            //! Note: all these values are in the NTT domain.
            lazy_muladd(evk_ct2[0].begin(), evk.beta_at(j).cptr_at(i), ct2.cbegin(), degree);
            lazy_muladd(evk_ct2[1].begin(), evk.alpha_at(j).cptr_at(i), ct2.cbegin(), degree);
        }
        // Lazy reduction, ct[l] + evk.alpha * ct[l] for l \in {0, 1}
        for (int l : {0, 1}) {
            auto rop_ptr = rop[l]->ptr_at(i);
            auto ct_ptr  = ct[l]->cptr_at(i);
            auto src_ptr = evk_ct2[l].cbegin();
            for (size_t d = 0; d < degree; ++d) {
                gT tmp = (*src_ptr++) + (*ct_ptr++);
                yell::ops::barret_reduction(&tmp, i);
                *rop_ptr++ = (T)tmp;
            }
        }
    }
}

// void Evaluator::Impl::Mul::rescale_ntt(
//   std::array<context::poly_t *, 2> rop,
//   std::array<context::poly_t *, 2> ct,
//   BaseConverter const& bconv) const
// {
//   constexpr size_t degree= context::degree;
//   long last_prime_index = ct[0]->moduli_count() - 1;
//   assert(last_prime_index >= 1);
//   for (int i : {0, 1}) {
//     assert(rop[i] && ct[i]);
//     assert(rop[i]->moduli_count() + 1 == ct[i]->moduli_count());
//     //! only convert the last moduli to power-basis
//     yell::ntt<degree>::backward(ct[i]->ptr_at(last_prime_index), last_prime_index);
//     bconv.rescale(rop[i], *(ct[i]));
//   }
// }

}  // namespace fHE
