#include "fHE/special_prime_chain.hpp"
#include <yell/ops.hpp>
#include "fHE/context.hpp"
namespace fHE {
SpecialPrimeChain::SpecialPrimeChain(size_t max_n_nrml_primes, size_t n_spcl_primes)
    : n_nrml_primes_(max_n_nrml_primes), n_spcl_primes_(n_spcl_primes) {
    yell::ops::mulmod mulmod;
    punc_prod_of_spcl_primes_.resize(max_n_nrml_primes);
    inv_prod_of_spcl_primes_.resize(max_n_nrml_primes);
    prod_of_spcl_primes_.resize(max_n_nrml_primes);
    for (size_t i = 0; i < max_n_nrml_primes; ++i) {
        punc_prod_of_spcl_primes_[i].resize(n_spcl_primes);
        for (size_t j = 0; j < n_spcl_primes; ++j) {
            T prod {1};
            for (size_t k = 0; k < n_spcl_primes; ++k) {
                if (j == k) continue;  // puncture
                size_t spcl_prime_idx = max_n_nrml_primes + k;
                mulmod.compute(prod, yell::params::P[spcl_prime_idx], i);
            }
            /// NOTE(juhou): we store the negated version so that we can use addmod instead of submod.
            punc_prod_of_spcl_primes_[i][j] = yell::params::P[i] - prod;
        }

        T prod_of_spcl_primes = yell::params::P[i] - punc_prod_of_spcl_primes_[i][0];
        T sp_0_idx            = max_n_nrml_primes;
        mulmod.compute(prod_of_spcl_primes, yell::params::P[sp_0_idx], i);
        prod_of_spcl_primes_[i]     = prod_of_spcl_primes;
        inv_prod_of_spcl_primes_[i] = yell::math::inv_mod_prime(prod_of_spcl_primes, i);
    }

    inv_punc_prod_of_spcl_primes_.resize(n_spcl_primes);
    for (size_t k = 0; k < n_spcl_primes; ++k) {
        T prod {1};
        size_t k_idx = max_n_nrml_primes + k;
        for (size_t d = 0; d < n_spcl_primes; ++d) {
            if (d == k) continue;
            size_t spcl_prime_idx = max_n_nrml_primes + d;
            mulmod.compute(prod, yell::params::P[spcl_prime_idx], k_idx);
        }
        inv_punc_prod_of_spcl_primes_[k] = yell::math::inv_mod_prime(prod, k_idx);
    }
}
}  // namespace fHE
