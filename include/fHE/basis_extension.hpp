#pragma once
#include <vector>
#include <memory>
#include <cstdint>
#include <unordered_map>
#include <yell/params.hpp>
#include "fHE/context.hpp"

namespace fHE {

class PrimeBundles;
class SpecialPrimeChain;

class BasisExtension {
public:
    using T      = yell::params::value_type;
    using gT     = yell::params::gt_value_type;
    using poly_t = context::poly_t;

    typedef struct ModupParams {
        /// Mod-up
        size_t bundle_idx;
        size_t n_nrml_primes;
        size_t n_spcl_primes;
    } ModupParams;

    typedef struct ModdownParams {
        size_t n_nrml_primes;
        size_t n_spcl_primes;
    } ModdownParams;

    explicit BasisExtension(std::shared_ptr<PrimeBundles> bundles, std::shared_ptr<SpecialPrimeChain> chain);

    ~BasisExtension();

    void mod_up(poly_t* out, const poly_t& in, ModupParams const& parms) const;

    void mod_down(poly_t* out, const poly_t& in, ModdownParams const& parms) const;

    void mod_down_inplace(poly_t* in, ModdownParams const& parms) const;
private:
    void mod_up_single(T* dst, const poly_t& in, std::vector<size_t> const& in_prime_idx, size_t dst_prime_idx) const;

    class NormalPrimeChain {
    public:
        typedef struct PunctureParam {
            size_t puncture_idx;
            size_t spcl_prime_idx;
        } PunctureParam;

        explicit NormalPrimeChain(size_t chain_length, size_t sp_idx_start, size_t n_spcl_primes);

        size_t chain_length() const { return chain_length_; }

        T inv_puncture_product(size_t nrml_prime_idx) const { return inv_punc_prods_.at(nrml_prime_idx); }

        T punc_prod_over_special_prime(PunctureParam const& parm) const {
            return punc_prods_MOD_spcl_prime_.at(parm.spcl_prime_idx).at(parm.puncture_idx);
        }

    private:
        size_t chain_length_;
        size_t sp_idx_start_;
        std::vector<T> inv_punc_prods_;
        std::vector<std::vector<T>> punc_prods_MOD_spcl_prime_;
    };

    const int max_n_nrml_primes_;
    std::shared_ptr<PrimeBundles> bundles_;
    std::shared_ptr<SpecialPrimeChain> specialPrimeChain_;
    std::unordered_map<size_t, std::shared_ptr<NormalPrimeChain>> normalPrimehains_;
};

}  // namespace fHE
