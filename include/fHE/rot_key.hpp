#pragma once
#include <memory>
#include <unordered_map>
#include <vector>
#include "fHE/context.hpp"
namespace fHE {
struct SK;  // forward declaration

struct SK;
class KSwithKeys;
class PrimeBundles;
class SpecialPrimeChain;

using poly_t = context::poly_t;

class RotKey {
public:
    enum class Direction { LEFT, RIGHT };

    explicit RotKey(SK const& sk,
                    uint16_t offset,
                    Direction dir,
                    std::shared_ptr<PrimeBundles> bundles,
                    std::shared_ptr<SpecialPrimeChain> chain);

    const poly_t& alpha_at(size_t j) const;

    const poly_t& beta_at(size_t j) const;

    size_t galois() const;

    ~RotKey();

private:
    std::shared_ptr<KSwithKeys> key_;
};

// struct MixedRotKey {
// private:
//     static constexpr size_t degree = context::degree;
//     static constexpr size_t L      = context::nr_ctxt_moduli;
//     using poly_t                   = context::poly_t;
//     std::vector<poly_t> alpha;  // L polynomials, each L + 1 moduli
//     std::vector<poly_t> beta;   // L polynomials, each L + 1 moduli
//
// public:
//     enum class Direction { LEFT, RIGHT };
//
//     const size_t galois;
//
//     explicit MixedRotKey(SK const& sk, uint16_t offset, Direction dir);
//
//     MixedRotKey(MixedRotKey const& oth) = delete;
//
//     MixedRotKey(MixedRotKey&& oth) = delete;
//
//     MixedRotKey& operator=(MixedRotKey const& oth) = delete;
//
//     const poly_t* get_beta_at(size_t j) const { return &beta.at(j); }
//
//     const poly_t* get_alpha_at(size_t j) const { return &alpha.at(j); }
//
//     poly_t* get_beta_at(size_t j) { return &beta.at(j); }
//
//     poly_t* get_alpha_at(size_t j) { return &alpha.at(j); }
//
//     ~MixedRotKey();
// };

template <class T>
struct KeyFactory {};

// template <>
// struct KeyFactory<MixedRotKey> {
//     std::shared_ptr<MixedRotKey> create(const SK& sk,
//                                         uint16_t idx,
//                                         MixedRotKey::Direction dir,
//                                         std::shared_ptr<PrimeBundles> #<{(||)}>#,
//                                         std::shared_ptr<SpecialPrimeChain> #<{(||)}>#) {
//         return std::make_shared<MixedRotKey>(sk, idx, dir);
//     }
// };

template <>
struct KeyFactory<RotKey> {
    std::shared_ptr<RotKey> create(const SK& sk,
                                   uint16_t idx,
                                   RotKey::Direction dir,
                                   std::shared_ptr<PrimeBundles> bundles,
                                   std::shared_ptr<SpecialPrimeChain> chain) {
        return std::make_shared<RotKey>(sk, idx, dir, bundles, chain);
    }
};

template <class RotKeyType>
struct RotKeySet {
private:
    static constexpr size_t logn = yell::static_log2<context::degree / 2>::value;
    using rotkey_ptr             = std::shared_ptr<RotKeyType>;
    std::unordered_map<size_t, rotkey_ptr> rot_keys;

public:
    RotKeySet(SK const& sk, std::shared_ptr<PrimeBundles> bundles, std::shared_ptr<SpecialPrimeChain> chain) {
        KeyFactory<RotKeyType> factory;
        for (size_t i = 0; i < logn; ++i) {
            int idx  = 1 << i;
            auto obj = factory.create(sk, idx, RotKeyType::Direction::LEFT, bundles, chain);
            rot_keys.insert({idx, std::move(obj)});
        }
    }

    RotKeySet(RotKeySet const& oth) = delete;

    RotKeySet(RotKeySet&& oth) = delete;

    RotKeySet& operator=(RotKeySet const& oth) = delete;

    const RotKeyType* get(uint16_t offset) const {
        if (offset >= (1LL << logn)) return nullptr;
        auto kv = rot_keys.find(offset);
        if (kv == rot_keys.end()) return nullptr;
        return kv->second.get();
    }

    ~RotKeySet() {}
};

}  // namespace fHE

