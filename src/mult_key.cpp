#include "fHE/keys.hpp"

#include "fHE/mult_key.hpp"
#include "fHE/kswitchkeys.hpp"
#include "fHE/prime_bundles.hpp"
#include "fHE/special_prime_chain.hpp"
#include "fHE/yell.hpp"

namespace fHE {

MultKey::MultKey(SK const& sk,
                 std::shared_ptr<PrimeBundles> const bundles,
                 std::shared_ptr<SpecialPrimeChain> const chain)
{
    const size_t n_spcl_primes = bundles->bundle_size();
    const size_t n_nrml_primes = bundles->max_n_nrml_primes();
    /// compute sx^2
    context::poly_t sx_square(sk.sx * sk.sx);

    key_ = std::make_shared<KSwithKeys>(sx_square, sk.sx, 1UL, bundles, chain);
    sx_square.clear();
}

MultKey::~MultKey() { }

const poly_t& MultKey::alpha_at(size_t j) const { return key_->alpha_at(j); }

const poly_t& MultKey::beta_at(size_t j) const { return key_->beta_at(j); }

} // namespace fHE
