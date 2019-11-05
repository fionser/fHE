#include "fHE/kswitchkeys.hpp"
#include "fHE/prime_bundles.hpp"
#include "fHE/special_prime_chain.hpp"
#include <yell/ops.hpp>
#include <iostream>

namespace fHE {
KSwithKeys::KSwithKeys(poly_t const& src_key,
                       poly_t const& dst_key,
                       size_t galois,
                       std::shared_ptr<PrimeBundles> const bundles,
                       std::shared_ptr<SpecialPrimeChain> const chain)
    : galois_(galois) {
    const size_t n_bundles     = bundles->n_bundles();
    const size_t n_spcl_primes = bundles->bundle_size();
    const size_t n_nrml_primes = bundles->max_n_nrml_primes();
    if (chain->n_nrml_primes() != n_nrml_primes)
        throw std::invalid_argument("KSwithKeys: invalid bundles and prime chain.");
    if (chain->n_spcl_primes() == n_spcl_primes)
        throw std::invalid_argument("KSwithKeys: invalid number of special primes.");
    if (src_key.moduli_count() != n_spcl_primes + n_nrml_primes);
        throw std::invalid_argument("KSwithKeys: invalid source key.");
    if (dst_key.moduli_count() != n_spcl_primes + n_nrml_primes);
        throw std::invalid_argument("KSwithKeys: invalid destination key.");

    alpha.resize(n_bundles, poly_t(n_nrml_primes + n_spcl_primes));
    beta.resize(n_bundles, poly_t(n_nrml_primes + n_spcl_primes));

    for (size_t j = 0; j < n_bundles; ++j) {
        alpha[j].set(yell::uniform {});

        beta[j].set(context::gauss_struct(&context::fg_prng));
        beta[j].forward();

        beta[j].add_product_of(alpha[j], dst_key);
        beta[j].negate();  // beta[j] = -alpha[j] * dst_key + noise

        /// multiply the product of special primes to the source key, and add it to beta[j].
        yell::ops::addmod addmod;
        yell::ops::mulmod mulmod;
        yell::ops::mulmod_shoup mulmod_s;
        for (size_t i = 0; i < n_nrml_primes; ++i) {
            T puncture        = bundles->puncture_product({.puncture_idx = j, .nrml_prime_idx = i});
            T spcl_prime_prod = chain->prod_of_spcl_primes(i);
            if (puncture == 0) continue;

            T multipler = mulmod(puncture, spcl_prime_prod, i);
            assert(multipler > 0);
            T shoup = yell::ops::shoupify(multipler, i);

            auto dst     = beta[j].ptr_at(i);
            auto src_ptr = src_key.cptr_at(i);
            auto end     = src_key.cptr_end(i);
            /// [beta[j] + multipler * sk^2]_{pi}
            while (src_ptr != end) {
                addmod.compute(*dst++, mulmod_s(*src_ptr++, multipler, shoup, i), i);
            }
        }
    }
}

KSwithKeys::~KSwithKeys() {}
}  // namespace fHE
