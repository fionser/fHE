#pragma once

#include <deque>
#include "fHE/context.hpp"

namespace fHE {

class PrimeBundles;
class SpecialPrimeChain;

class KSwithKeys {
public:
    using poly_t = context::poly_t;
    using T      = yell::params::value_type;

    /// Create a switching key from the source key (src_sx) to the destination key (dst_key)
    explicit KSwithKeys(poly_t const& src_key,
                        poly_t const& dst_key,
                        size_t galois,
                        std::shared_ptr<PrimeBundles> const bundles,
                        std::shared_ptr<SpecialPrimeChain> const chain);

    ~KSwithKeys();

    KSwithKeys(KSwithKeys&& oth) : galois_(oth.galois_), alpha(std::move(oth.alpha)), beta(std::move(oth.beta)) {}

    KSwithKeys(const KSwithKeys& oth) : galois_(oth.galois_), alpha(oth.alpha), beta(oth.beta) {}

    KSwithKeys& operator=(const KSwithKeys& oth) = delete;

    const poly_t& alpha_at(size_t j) const { return alpha.at(j); }

    const poly_t& beta_at(size_t j) const { return beta.at(j); }

    const size_t& galois() const { return galois_; }

protected:
    size_t galois_;
    std::deque<poly_t> alpha;
    std::deque<poly_t> beta;
};

};  // namespace fHE
