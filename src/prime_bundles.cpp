#include "fHE/prime_bundles.hpp"
#include <cmath>
#include <numeric>
#include <yell/ops.hpp>
/* ceil(a / b) */
static inline size_t div_ceil(size_t a, size_t b) {
    assert(b > 0);
    return (a + b - 1) / b;
}

namespace fHE {

PrimeBundles::PrimeBundles(size_t max_n_nrml_primes, size_t bundle_size)
    : max_n_nrml_primes_(max_n_nrml_primes),
      //bundle_size_(std::max(1UL, static_cast<size_t>(std::floor(std::sqrt(max_n_nrml_primes * 1.))))),
      //bundle_size_(1), /// SEAL v3.2
      //bundle_size_(max_n_nrml_primes), /// HEAAN
      bundle_size_(bundle_size),
      n_bundles_(div_ceil(max_n_nrml_primes, bundle_size_)),
      bundled_prods_(std::vector<std::vector<T>>(n_bundles_, std::vector<T>(max_n_nrml_primes, 0))) {
    assert(max_n_nrml_primes > 0 && max_n_nrml_primes < yell::params::kMaxNbModuli);
    assert(bundle_size_ * n_bundles_ >= max_n_nrml_primes);
    assert(bundle_size_ >= 1 && n_bundles_ >= 1);

    for (size_t prime_idx = 0; prime_idx < max_n_nrml_primes; ++prime_idx) {
        for (size_t bundle_idx = 0; bundle_idx < n_bundles_; ++bundle_idx) {
            size_t prime_start = bundle_idx * bundle_size_;
            size_t prime_end   = std::min(prime_start + bundle_size_, max_n_nrml_primes);
            if (prime_start <= prime_idx && prime_idx < prime_end) {
                // this prime inside the current bundle, so the bundled product should be 0.
                continue;
            }

            T product {1};
            yell::ops::mulmod mulmod;
            for (size_t j = prime_start; j < prime_end; ++j) {
                mulmod.compute(product, yell::params::P[j], prime_idx);
            }
            bundled_prods_.at(bundle_idx).at(prime_idx) = product;
        }
    }
}

std::vector<size_t> PrimeBundles::prime_indices_in_bundle(size_t bundle_idx, size_t n_nrml_primes) const {
    assert(n_nrml_primes <= max_n_nrml_primes());
    if (n_nrml_primes == 0) n_nrml_primes = max_n_nrml_primes();
    std::vector<size_t> indices;
    if (bundle_idx >= n_bundles()) return indices;
    size_t bsze     = bundle_size();
    size_t n_primes = std::min((bundle_idx + 1) * bsze, n_nrml_primes) - bundle_idx * bsze;
    indices.resize(n_primes);
    std::iota(indices.begin(), indices.end(), bundle_idx * bsze);
    return indices;
}

PrimeBundles::T PrimeBundles::puncture_product(PunctureParam const& parm) const {
    assert(parm.puncture_idx < bundled_prods_.size());

    T prod {1};
    const size_t nrml_prime_idx = parm.nrml_prime_idx;
    yell::ops::mulmod mulmod;
    for (size_t bundle_idx = 0; bundle_idx < n_bundles_; ++bundle_idx) {
        if (bundle_idx == parm.puncture_idx) continue;  // puncture
        T bundled_prod = bundled_product({.bundle_idx = bundle_idx, .nrml_prime_idx = nrml_prime_idx});
        mulmod.compute(prod, bundled_prod, nrml_prime_idx);
    }
    return prod;
}

}  // namespace fHE
