#include "fHE/rot_key.hpp"
#include "fHE/keys.hpp"
#include "fHE/kswitchkeys.hpp"
#include "fHE/prime_bundles.hpp"
#include "fHE/special_prime_chain.hpp"
#include "fHE/yell.hpp"

#include "fHE/encoder.hpp"

namespace fHE {

static size_t rotation_offset_galois(int offset) {
    constexpr size_t nr_slots  = context::degree >> 1u;
    constexpr size_t m         = context::degree << 1u;
    constexpr size_t generator = Encoder::generator;
    bool sign                  = offset > 0;
    offset                     = std::abs(offset) % nr_slots;
    if (sign) offset = nr_slots - offset;
    size_t galois {1};
    while (offset--) {  // g^{k} mod m
        galois *= generator;
        galois &= (m - 1);
    }
    return galois;
}

RotKey::RotKey(SK const &sk,
               uint16_t offset,
               Direction dir,
               std::shared_ptr<PrimeBundles> bundles,
               std::shared_ptr<SpecialPrimeChain> chain) {

    const size_t n_spcl_primes = bundles->bundle_size();
    const size_t n_nrml_primes = bundles->max_n_nrml_primes();
    assert(sk.sx.moduli_count() == n_spcl_primes + n_nrml_primes);
    size_t galois = rotation_offset_galois(dir == Direction::RIGHT ? offset : -(int) offset);

    context::poly_t galois_sx(sk.sx);
    apply_galois(&galois_sx, galois);
    key_ = std::make_shared<KSwithKeys>(galois_sx, sk.sx, galois, bundles, chain);
}

RotKey::~RotKey() {}

const poly_t& RotKey::alpha_at(size_t j) const { return key_->alpha_at(j); }

const poly_t& RotKey::beta_at(size_t j) const { return key_->beta_at(j); }

size_t RotKey::galois() const { return key_->galois(); }

#if 0
MixedRotKey::MixedRotKey(SK const &sk, uint16_t offset, Direction dir)
    : galois(rotation_offset_galois(dir == Direction::RIGHT ? offset : -(int)offset)),
      beta(std::vector<poly_t>(L, poly_t(L + 1))),
      alpha(std::vector<poly_t>(L, poly_t(L + 1, yell::uniform {}))) {
    using T                 = yell::params::value_type;
    constexpr size_t degree = context::degree;
    constexpr size_t bytes  = degree * sizeof(T);

    context::poly_t ext_sx(L + 1);
    for (size_t j = 0; j < L + 1; ++j) {
        std::memcpy(ext_sx.ptr_at(j), sk.sx.cptr_at(j), bytes);
    }

    auto rotated_sx(ext_sx);
    apply_galois(&rotated_sx, galois);

    // beta_j = -sx * alpha_j + [qK]_qj * rotated(sx) + e
    yell::ops::muladd muladd;
    T qK = yell::params::P[context::index_sp_prime(0)];
    for (size_t j = 0; j < L; ++j) {
        context::poly_t *beta_        = get_beta_at(j);
        const context::poly_t *alpha_ = get_alpha_at(j);
        beta_->set(context::gauss_struct(&context::fg_prng));
        beta_->forward();

        beta_->add_product_of(*alpha_, ext_sx);
        beta_->negate();  // -e - alpha * sx

        std::transform(beta_->cptr_at(j), beta_->cptr_end(j), rotated_sx.cptr_at(j), beta_->ptr_at(j),
                       [muladd, j, qK](T const &beta, T const &s) -> T {
                           // [beta]_qj += qK * [rotated_sx]_qj
                           return muladd(beta, s, qK, j);
                       });
    }
    rotated_sx.clear();
    ext_sx.clear();
}

MixedRotKey::~MixedRotKey() {}
#endif
}  // namespace fHE
