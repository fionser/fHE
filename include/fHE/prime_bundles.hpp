#pragma once
#include <vector>
#include <cstdint>
#include <yell/params.hpp>

namespace fHE {

class PrimeBundles {
public:
    using T = yell::params::value_type;

    typedef struct PunctureParam {
        size_t puncture_idx;
        size_t nrml_prime_idx;
    } PunctureParam;

    typedef struct BunldedParam {
        size_t bundle_idx;
        size_t nrml_prime_idx;
    } BunldedParam;

    explicit PrimeBundles(size_t max_n_nrml_primes, size_t n_spcl_primes);

    std::vector<size_t> prime_indices_in_bundle(size_t bundle_idx, size_t n_nrml_primes = 0) const;

    /// To get a bundled product (i.e., the product of all primes in the specified bundle),
    /// over the specified normal prime.
    inline T bundled_product(BunldedParam const& parm) const {
        assert(parm.bundle_idx < n_bundles_);
        assert(parm.nrml_prime_idx < max_n_nrml_primes_);
        return bundled_prods_.at(parm.bundle_idx).at(parm.nrml_prime_idx);
    }

    T puncture_product(PunctureParam const& parm) const;

    inline size_t max_n_nrml_primes() const { return max_n_nrml_primes_; }

    inline size_t bundle_size() const { return bundle_size_; }

    inline size_t n_bundles() const { return n_bundles_; }

private:
    size_t max_n_nrml_primes_;
    size_t bundle_size_;  // ceil(max_n_nrml_primes / n_bundles)
    size_t n_bundles_;    // floor(sqrt(max_n_nrml_primes)))
    std::vector<std::vector<T>> bundled_prods_;
};
} // namespace fHE
