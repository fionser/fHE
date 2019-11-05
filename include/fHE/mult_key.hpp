#pragma once
#include <memory>
#include "fHE/context.hpp"

namespace fHE {

struct SK;
class KSwithKeys;
class PrimeBundles;
class SpecialPrimeChain;
using poly_t = context::poly_t;

class MultKey {
public:
    explicit MultKey(SK const& sk,
                     std::shared_ptr<PrimeBundles> const bundles,
                     std::shared_ptr<SpecialPrimeChain> const chain);

    const poly_t& alpha_at(size_t j) const;

    const poly_t& beta_at(size_t j) const;

    constexpr size_t galois() const { return 1UL; }

    ~MultKey();

private:
    std::shared_ptr<KSwithKeys> key_;
};
}  // namespace fHE

