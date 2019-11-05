#include "fHE/basis_extension.hpp"
#include "fHE/prime_bundles.hpp"
#include "fHE/special_prime_chain.hpp"

/* ceil(a / b) */
static inline size_t div_ceil(size_t a, size_t b) {
    assert(b > 0);
    return (a + b - 1) / b;
}

namespace fHE {
BasisExtension::BasisExtension(std::shared_ptr<PrimeBundles> bundles, std::shared_ptr<SpecialPrimeChain> chain)
    : max_n_nrml_primes_(bundles->max_n_nrml_primes()), bundles_(bundles), specialPrimeChain_(chain) {
    const size_t n_nrml_primes = bundles->max_n_nrml_primes();
    const size_t n_spcl_primes = bundles->bundle_size();
    for (size_t chain_length = 1; chain_length <= n_nrml_primes; ++chain_length) {
        auto normal = std::make_shared<NormalPrimeChain>(chain_length, max_n_nrml_primes_, n_spcl_primes);
        normalPrimehains_.insert({chain_length, std::move(normal)});
    }
}

BasisExtension::~BasisExtension() {}

void BasisExtension::mod_up(poly_t* out, const poly_t& in, ModupParams const& parms) const {
    if (!out) return;

    std::vector<size_t> cur_indices = bundles_->prime_indices_in_bundle(parms.bundle_idx);
    const size_t cur_n_bundles      = div_ceil(parms.n_nrml_primes, bundles_->bundle_size());
    if (parms.bundle_idx >= cur_n_bundles)
        throw std::invalid_argument("BasisExtension::mod_up: parms.bundle_idx too large.");
    if (parms.n_nrml_primes > bundles_->max_n_nrml_primes())
        throw std::invalid_argument("BasisExtension::mod_up: parms.n_nrml_primes too large.");
    if (in.moduli_count() > cur_indices.size())
        throw std::invalid_argument("BasisExtension::mod_up: in.moduli_count() too large.");
    if (out->moduli_count() != parms.n_nrml_primes + parms.n_spcl_primes)
        throw std::invalid_argument("BasisExtension::mod_up: invalid out->moduli_count().");

    /// The number of current normal primes might be smaller than the bundle size.
    /// Drop the lastest few primes.
    while (cur_indices.size() > in.moduli_count()) {
        cur_indices.pop_back();
    }

    for (size_t bundle_idx = 0; bundle_idx < cur_n_bundles; ++bundle_idx) {
        if (bundle_idx == parms.bundle_idx) {
            const size_t nbytes = sizeof(in(0, 0)) * in.degree;
            for (size_t _i = 0; _i < cur_indices.size(); ++_i) {
                std::memcpy(out->ptr_at(cur_indices[_i]), in.cptr_at(_i), nbytes);
            }
        } else {
            std::vector<size_t> dst_indices = bundles_->prime_indices_in_bundle(bundle_idx);
            for (size_t dst_prime_index : dst_indices) {
                if (dst_prime_index >= parms.n_nrml_primes) continue;
                mod_up_single(out->ptr_at(dst_prime_index), in, cur_indices, dst_prime_index);
            }
        }
    }

    for (size_t k = 0; k < parms.n_spcl_primes; ++k) {
        mod_up_single(out->ptr_at(parms.n_nrml_primes + k), in, cur_indices, max_n_nrml_primes_ + k);
    }
}

void BasisExtension::mod_down_inplace(poly_t* in, ModdownParams const& parms) const {
    /// The special prime part of [in] is in the power-basis.
    /// But, the normal prime part of [in] is in the NTT-basis.
    if (!in || (in->moduli_count() != parms.n_spcl_primes + parms.n_nrml_primes)) {
        throw std::invalid_argument("BasisExtension::mod_down_inplace: invalid input parameters");
    }

    const size_t degree = in->degree;
    yell::ops::mulmod_shoup mulmod_s;
    std::vector<gT> lazy(degree);
    std::vector<T> ntt_temp(degree);
    for (size_t nrml_idx = 0; nrml_idx < parms.n_nrml_primes; ++nrml_idx) {
        std::memset(lazy.data(), 0, sizeof(gT) * degree);
        /// Step 1: Conv_{Specail Primes -> Normal Primes}
        ///        [c']_{q_i} = \sum_j { { [c]_{p_j} * \hat{p}_j^{-1}  mod p_j } * (-\hat{p}_j) mod q_i }
        ///        Note: q_i is normal prime, and p_j is special prime.
        for (size_t spcl_idx = 0; spcl_idx < parms.n_spcl_primes; ++spcl_idx) {
            const size_t special_prime   = max_n_nrml_primes_ + spcl_idx;
            auto in_ptr                  = in->cptr_at(parms.n_nrml_primes + spcl_idx);
            auto dst_ptr                 = lazy.begin();
            const T inv_punch_prod       = specialPrimeChain_->inv_punc_prod_of_spcl_primes(spcl_idx);
            const T inv_punch_prod_shoup = yell::ops::shoupify(inv_punch_prod, special_prime);
            /// Opt: we use the negative puncture product, so that we use addition which is faster than subtraction.
            const gT neg_punch_prod = specialPrimeChain_->neg_punc_prod_of_spcl_primes(
                {.nrml_prime_idx = nrml_idx, .puncture_idx = spcl_idx});

            while (dst_ptr != lazy.end()) {
                T v = mulmod_s(*in_ptr++, inv_punch_prod, inv_punch_prod_shoup, special_prime);
                (*dst_ptr++) += v * neg_punch_prod;
            }
        }

        /// Step 2: lazy reduction
        std::transform(lazy.cbegin(), lazy.cend(), ntt_temp.data(), [nrml_idx](gT v) -> T {
            yell::ops::barret_reduction(&v, nrml_idx);
            return (T)v;
        });

        /// Step 3: forward the normal prime part
        yell::ntt<context::degree>::forward(ntt_temp.data(), nrml_idx);

        /// Step 4: [P^{-1}]_{q_i} * ([c]_{q_i} + [c']_{q_i}) mod q_i
        const T P       = specialPrimeChain_->inv_prod_of_spcl_primes(nrml_idx);
        const T P_shoup = yell::ops::shoupify(P, nrml_idx);
        std::transform(ntt_temp.cbegin(), ntt_temp.cend(), in->cptr_at(nrml_idx), in->ptr_at(nrml_idx),
                       [mulmod_s, P, P_shoup, nrml_idx](const T a, const T b) -> T {
                           return mulmod_s(a + b, P, P_shoup, nrml_idx);
                       });
    }
}

void BasisExtension::mod_down(poly_t* out, const poly_t& in, ModdownParams const& parms) const {
    /// The special prime part of [in] is in the power-basis.
    /// But, the normal prime part of [in] is in the NTT-basis.
    /// By using this optimization, we can save some NTTs & invNTTs.
    if (!out) return;
    if ((out->moduli_count() != parms.n_nrml_primes) ||
        (in.moduli_count() != parms.n_spcl_primes + parms.n_nrml_primes)) {
        throw std::invalid_argument("BasisExtension::mod_down: invalid input parameters");
    }

    yell::ops::mulmod mulmod;
    for (size_t nrml_idx = 0; nrml_idx < parms.n_nrml_primes; ++nrml_idx) {
        std::vector<gT> lazy(in.degree, 0);
        for (size_t spcl_idx = 0; spcl_idx < parms.n_spcl_primes; ++spcl_idx) {
            size_t special_prime = max_n_nrml_primes_ + spcl_idx;
            auto in_ptr          = in.cptr_at(parms.n_nrml_primes + spcl_idx);
            auto dst_ptr         = lazy.data();
            T inv_punch_prod     = specialPrimeChain_->inv_punc_prod_of_spcl_primes(spcl_idx);
            /// Opt: we use the negative puncture product, so that we use addition which is faster than subtraction.
            T neg_punch_prod = specialPrimeChain_->neg_punc_prod_of_spcl_primes(
                {.nrml_prime_idx = nrml_idx, .puncture_idx = spcl_idx});
            for (size_t d = 0; d < in.degree; ++d, ++dst_ptr) {
                T v = mulmod(*in_ptr++, inv_punch_prod, special_prime);
                (*dst_ptr) += (gT)v * neg_punch_prod;
            }
        }

        std::transform(lazy.cbegin(), lazy.cend(), out->ptr_at(nrml_idx), [nrml_idx](gT v) -> T {
            yell::ops::barret_reduction(&v, nrml_idx);
            return (T)v;
        });
        yell::ntt<context::degree>::forward(out->ptr_at(nrml_idx), nrml_idx);
    }

    /// We have used the negative puncture product, so that we need addition instead of subtraction here.
    for (size_t nrml_idx = 0; nrml_idx < parms.n_nrml_primes; ++nrml_idx) {
        yell::ops::mulmod_shoup mulmod_s;
        T P     = specialPrimeChain_->inv_prod_of_spcl_primes(nrml_idx);
        T shoup = yell::ops::shoupify(P, nrml_idx);
        std::transform(in.cptr_at(nrml_idx), in.cptr_end(nrml_idx), out->cptr_at(nrml_idx), out->ptr_at(nrml_idx),
                       [mulmod_s, P, shoup, nrml_idx](const T a, const T b) -> T {
                           // addition with lazy reduction
                           return mulmod_s(a + b, P, shoup, nrml_idx);
                       });
    }
}

void BasisExtension::mod_up_single(T* dst,
                                   poly_t const& in,
                                   std::vector<size_t> const& in_prime_idx,
                                   size_t dst_prime_idx) const {
    if (!dst) return;
    assert(in.moduli_count() == in_prime_idx.size());
    assert(dst_prime_idx < yell::params::kMaxNbModuli);
    for (size_t idx : in_prime_idx) {
        assert(idx < bundles_->max_n_nrml_primes());
        assert(idx != dst_prime_idx);
    }

    std::vector<T> inv_punch_prod;
    std::vector<T> punch_prod;
    yell::ops::mulmod mulmod;
    for (size_t punch_idx : in_prime_idx) {
        T inv_prod {1}, prod {1};
        for (size_t prime_idx : in_prime_idx) {
            if (punch_idx == prime_idx) continue;
            mulmod.compute(prod, yell::params::P[prime_idx], dst_prime_idx);
            mulmod.compute(inv_prod, yell::params::P[prime_idx], punch_idx);
        }
        punch_prod.push_back(prod);
        inv_punch_prod.push_back(yell::math::inv_mod_prime(inv_prod, punch_idx));
    }

    yell::ops::mulmod_shoup mulmod_s;
    std::vector<gT> accum(in.degree, 0);  // use lazy reduction
    for (size_t i = 0; i < in.moduli_count(); ++i) {
        size_t prime_idx = in_prime_idx[i];
        T shoup          = yell::ops::shoupify(inv_punch_prod.at(i), prime_idx);

        const T* src_ptr = in.cptr_at(i);
        auto dst_ptr = accum.begin();
        while (dst_ptr != accum.end()) {
            T tmp = mulmod_s(*src_ptr++, inv_punch_prod.at(i), shoup, prime_idx);
            (*dst_ptr++) += (gT) tmp * punch_prod[i];
        }
    }

    std::transform(accum.cbegin(), accum.cend(), dst, [dst_prime_idx](gT v) -> T {
        yell::ops::barret_reduction(&v, dst_prime_idx);
        return (T)v;
    });
}

BasisExtension::NormalPrimeChain::NormalPrimeChain(size_t chain_length, size_t sp_idx_start, size_t n_spcl_primes)
    : chain_length_(chain_length), sp_idx_start_(sp_idx_start) {
    assert(chain_length > 0);
    assert(n_spcl_primes > 0);
    yell::ops::mulmod mulmod;

    inv_punc_prods_.resize(chain_length);
    for (size_t i = 0; i < chain_length; ++i) {
        T prod {1};
        for (size_t j = 0; j < chain_length; ++j) {
            if (j == i) continue;
            mulmod.compute(prod, yell::params::P[j], i);
        }
        inv_punc_prods_[i] = yell::math::inv_mod_prime(prod, i);
    }

    punc_prods_MOD_spcl_prime_.resize(n_spcl_primes);
    for (size_t k = 0; k < n_spcl_primes; ++k) {
        size_t spcl_idx = sp_idx_start_ + k;
        punc_prods_MOD_spcl_prime_[k].resize(chain_length);
        for (size_t i = 0; i < chain_length; ++i) {
            T prod {1};
            for (size_t j = 0; j < chain_length; ++j) {
                if (j == i) continue;
                mulmod.compute(prod, yell::params::P[j], spcl_idx);
            }

            punc_prods_MOD_spcl_prime_.at(k).at(i) = prod;
        }
    }
}

}  // namespace fHE
