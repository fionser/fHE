#pragma once
#include <cstdint>
#include <vector>
#include <yell/params.hpp>

namespace fHE {
class SpecialPrimeChain {
public:
    using T = yell::params::value_type;

    typedef struct PunctureParam {
        size_t puncture_idx;
        size_t nrml_prime_idx;
    } PunctureParam;

    explicit SpecialPrimeChain(size_t max_n_nrml_primes, size_t n_spcl_primes);

    inline size_t n_nrml_primes() const { return n_nrml_primes_; }

    inline size_t n_spcl_primes() const { return n_spcl_primes_; }

    T inv_punc_prod_of_spcl_primes(size_t spcl_prime_idx) const {
        return inv_punc_prod_of_spcl_primes_.at(spcl_prime_idx);
    }

    T inv_prod_of_spcl_primes(size_t nrml_prime_idx) const { return inv_prod_of_spcl_primes_.at(nrml_prime_idx); }

    T prod_of_spcl_primes(size_t nrml_prime_idx) const { return prod_of_spcl_primes_.at(nrml_prime_idx); }

    T neg_punc_prod_of_spcl_primes(PunctureParam const& parm) const {
        return punc_prod_of_spcl_primes_.at(parm.nrml_prime_idx).at(parm.puncture_idx);
    }

private:
    size_t n_nrml_primes_;
    size_t n_spcl_primes_;
    std::vector<std::vector<T>> punc_prod_of_spcl_primes_;
    std::vector<T> inv_punc_prod_of_spcl_primes_;
    std::vector<T> prod_of_spcl_primes_;
    std::vector<T> inv_prod_of_spcl_primes_;
};
}  // namespace fHE
