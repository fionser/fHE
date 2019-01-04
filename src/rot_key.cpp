#include "fHE/rot_key.hpp"
#include "fHE/keys.hpp"
#include "fHE/yell.hpp"
namespace fHE {
RotKey::RotKey(SK const &sk, uint16_t offset, Direction dir)
    : galois(rotation_offset_galois(dir == Direction::RIGHT ? offset : - (int) offset)),
       beta(L + K),
      alpha(L + K, yell::uniform{})
{
  using  T = yell::params::value_type;
  using gT = yell::params::gt_value_type;
  //! At the end, we compute beta = P * s' - e - alpha * s
  //! s' is the rotated secret key, and P is the product of all special moduli.
  auto rotated_sx(sk.sx);
  apply_galois(&rotated_sx, galois);
  lift_normal_moduli(&rotated_sx);
  //! Reuse beta as the 'e'.
  //! e + alpha * s (with special moduli) using lazy reduction
  beta.set(context::gauss_struct(&context::fg_prng));
  beta.forward_lazy();
  for (size_t cm = 0; cm < L + K; ++cm) {
    auto bptr = beta.ptr_at(cm);
    auto aptr = alpha.cptr_at(cm);
    auto sptr = cm < L ? sk.sx.cptr_at(cm) : sk.ext_sx.cptr_at(cm - L);
    for (size_t d = 0; d < degree; ++d) {
      gT tmp = (gT) (*aptr++) * (*sptr++) + *bptr;
      yell::ops::barret_reduction(&tmp, cm);
      *bptr++ = (T) tmp;
    }
  }
  beta.negate();
  //! Finally, adding P * s' to beta
  yell::ops::addmod addmod;
  for (size_t cm = 0; cm < L; ++cm) {
    auto op0 = beta.ptr_at(cm);
    auto op1 = rotated_sx.cptr_at(cm);
    for (size_t d = 0; d < degree; ++d)
      addmod.compute(*op0++, *op1++, cm);
  }

  rotated_sx.clear();
}

RotKey::~RotKey() {}

size_t RotKey::rotation_offset_galois(int offset) const
{
  constexpr size_t nr_slots = context::degree >> 1u;
  constexpr size_t m = context::degree << 1u;
  constexpr size_t generator = 3;
  bool sign = offset > 0;
  offset = std::abs(offset) % nr_slots;
  if (sign)
    offset = nr_slots - offset;
  size_t galois{1};
  while(offset--) { // g^{k} mod m
    galois *= generator;
    galois &= (m - 1);
  }
  return galois;
}

void RotKey::lift_normal_moduli(context::poly_t *op) const
{
  if (!op) return;
  assert(op->moduli_count() == L);
  yell::ops::mulmod mulmod;
  std::vector<size_t> sp_primes(K);
  std::iota(sp_primes.begin(), sp_primes.end(),
            context::index_sp_prime(0));
  for (size_t cm = 0; cm < L; ++cm) {
    /* multiply P to each normal primes */
    auto P_mod_qj = product_of(sp_primes.begin(), sp_primes.end(), cm);
    auto dst = op->ptr_at(cm);
    for (size_t d = 0; d < degree; ++d)
      mulmod.compute(*dst++, P_mod_qj, cm);
  }
}
}// namespace fHE
