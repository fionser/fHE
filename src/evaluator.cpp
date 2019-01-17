#include "fHE/evaluator.hpp"
#include "fHE/cipher.hpp"
#include "fHE/rot_key.hpp"
#include "fHE/mult_key.hpp"
#include "fHE/base_converter.hpp"
#include "fHE/encoder.hpp"
#include "fHE/yell.hpp"
#include <array>

namespace fHE {

struct Evaluator::Impl {
  BaseConverter bconv;

  inline bool can_add_with(Cipher const& c0, Cipher const& c1) const {
    // TODO check scale
    if (c0.moduli_count() != c1.moduli_count())
      return false;
    return true;
  }

  inline bool can_mul_with(Cipher const& c0, Cipher const& c1) const {
    if (c0.moduli_count() != c1.moduli_count())
      return false;
    if (!(c0.canonical() && c1.canonical()))
      return false;
    return true;
  }

  bool add_inplace(Cipher *rop, Cipher const& c1) const;

  bool multiply_inplace(Cipher *rop, Cipher const& c1) const;

  bool relinerize_inplace(Cipher *rop, MultKey const& mkey) const;

  bool rotate_slots_inplace(Cipher*rop, RotKey const& rkey) const;

  bool rotate_slots_inplace_non_ntt(Cipher*rop, RotKey const& rkey) const;

private:
  struct Mul {
    void cross_mult(context::poly_t ans[3],
                    Cipher const& c0,
                    Cipher const& c1) const;

    void reLinearize(std::array<context::poly_t *, 2> &rop,
                     std::array<const context::poly_t *, 3> ct,
                     MultKey const& evk) const;

    void rescale_ntt(std::array<context::poly_t *, 2> rop,
                     std::array<context::poly_t *, 2> ct,
                     BaseConverter const& bconv) const;
  };

  void apply_rotation_key(std::array<context::poly_t *, 2> rop,
                          std::array<context::poly_t *, 2> d,
                          RotKey const &rotk) const;

  void apply_rotation_key_non_ntt(std::array<context::poly_t *, 2> rop,
                                  std::array<context::poly_t *, 2> d,
                                  RotKey const &rotk) const;
};

Evaluator::Evaluator(std::shared_ptr<Encoder> encoder) 
  : impl_(std::make_shared<Impl>()),
    encoder(encoder) {}

Evaluator::~Evaluator() {}

bool Evaluator::add(Cipher *rop, Cipher const& c1) const
{
  if (!impl_) return false;
  return impl_->add_inplace(rop, c1);
}

bool Evaluator::multiply(Cipher *rop, Cipher const& c1) const
{
  if (!impl_) return false;
  return impl_->multiply_inplace(rop, c1);
}

bool Evaluator::multiply_plain(Cipher *rop, 
                               context::poly_t const& p,
                               const double plain_scale) const
{
  if (!rop) return false;
  size_t nm = p.moduli_count();
  if (rop->moduli_count() != nm)
    return false;
  yell::ops::mulmod mulmod;
  for (size_t cm = 0; cm < rop->moduli_count(); ++cm) {
    auto op0 = rop->bx.ptr_at(cm);
    auto op1 = rop->ax.ptr_at(cm);
    auto op2 = p.cptr_at(cm);
    for (size_t d = 0; d < context::degree; ++d) {
      auto _op2 = *op2++;
      mulmod.compute(*op0++, _op2, cm);
      mulmod.compute(*op1++, _op2, cm);
    }
  }

  if (plain_scale > 0.)
    rop->scale(rop->scale() * plain_scale);
  return true;
}

bool Evaluator::relinearize(Cipher *rop, MultKey const& mkey) const
{
  if (!impl_) return false;
  return impl_->relinerize_inplace(rop, mkey);
}

bool Evaluator::multiply_and_relin(Cipher *rop, Cipher const& c1, MultKey const& mkey) const
{
  if (!impl_) return false;
  if (impl_->multiply_inplace(rop, c1))
    return impl_->relinerize_inplace(rop, mkey);
  else
    return false;
}

bool Evaluator::rotate_slots(Cipher*rop, RotKey const& rkey) const
{
  if (!impl_ || !rop) return false;
  if (!rop->canonical()) {
    std::cerr << "Can only rotate canonical ciphertext, call relinearize first." << std::endl;
    return false;
  }
  return impl_->rotate_slots_inplace(rop, rkey);
}

bool Evaluator::rotate_slots_non_ntt(Cipher*rop, RotKey const& rkey) const
{
  if (!impl_ || !rop) return false;
  if (!rop->canonical()) {
    std::cerr << "Can only rotate canonical ciphertext, call relinearize first." << std::endl;
    return false;
  }
  return impl_->rotate_slots_inplace_non_ntt(rop, rkey);
}

bool Evaluator::rotate_slots(Cipher *rop, int offset, RotKeySet const& rkeys) const
{
  if (!impl_ || !rop) return false;
  if (!rop->canonical()) {
    std::cerr << "Can only rotate canonical ciphertext, call relinearize first." << std::endl;
    return false;
  }
  constexpr size_t nr_slots = context::degree >> 1u;
  offset = std::abs(offset) % nr_slots;
  int steps = (int) (std::log2((double) offset) + 1);
  for (int i = 0; i < steps; ++i) {
    if (offset & (1 << i)) {
      RotKey const* rotkey = rkeys.get(1 << i);
      if (!rotkey) {
        std::cerr << "No such rotation key for offset " << (1 << i) << "\n";
        return false;
      }
      rotate_slots(rop, *rotkey);
    }
  }
  return true;
}

bool Evaluator::rotate_slots_non_ntt(Cipher *rop, int offset, RotKeySet const& rkeys) const
{
  if (!impl_ || !rop) return false;
  if (!rop->canonical()) {
    std::cerr << "Can only rotate canonical ciphertext, call relinearize first." << std::endl;
    return false;
  }
  constexpr size_t nr_slots = context::degree >> 1u;
  offset = std::abs(offset) % nr_slots;
  const int steps = (int) (std::log2((double) offset) + 1);
  for (int i = 0; i < steps; ++i) {
    if (offset & (1 << i)) {
      RotKey const* rotkey = rkeys.get(1 << i);
      if (!rotkey) {
        std::cerr << "No such rotation key for offset " << (1 << i) << "\n";
        return false;
      }
      rotate_slots_non_ntt(rop, *rotkey);
    }
  }
  return true;
}

bool Evaluator::replicate(Cipher *rop, size_t pos, RotKeySet const &rkeys) const
{
  if (!impl_ || !rop || !encoder) return false;
  if (!rop->canonical()) {
    std::cerr << "Can only replicate canonical ciphertext, call relinearize first." << std::endl;
    return false;
  }
  constexpr size_t nr_slots = context::degree >> 1;
  constexpr size_t logn = yell::static_log2<nr_slots>::value;
  //! replicate uses one masking and logN rotations and additions.
  pos %= (nr_slots);
  std::vector<double> mask(nr_slots, 0.);
  mask[pos] = 1.;

  double scale = context::encoder_scale / 8.;
  context::poly_t mask_poly(rop->moduli_count());
  encoder->encode(&mask_poly, mask, scale);
  mask_poly.forward_lazy();
  multiply_plain(rop, mask_poly, scale);
  //! rotate above the power-basis is faster
  rop->bx.backward();
  rop->ax.backward();

  for (size_t i = 0; i < logn; ++i) {
    fHE::Cipher tmp(*rop);
    rotate_slots_non_ntt(&tmp, 1 << i, rkeys);
    add(rop, tmp); // add is also ok on the power-basis.
  }

  rop->bx.forward();
  rop->ax.forward();
  return true;
}

// Implementation details go from here
bool Evaluator::Impl::add_inplace(
  Cipher *rop,
  Cipher const& c1) const
{
  if (!rop) return false;
  if (!can_add_with(*rop, c1)) {
    std::cerr << "Can not add_inplace." << std::endl;
    return false;
  }

  if (!rop->canonical() && c1.canonical()) {
    std::cerr << "Can not add_inplace" << std::endl;
    return false;
  }
  rop->bx += c1.bx;
  rop->ax += c1.ax;

  if (!c1.canonical()) {
    if (rop->canonical())
      rop->mul_aux = std::make_shared<context::poly_t>(*c1.mul_aux);
    else
      (*rop->mul_aux) += (*c1.mul_aux);
  } 

  if (rop->scale() == 1.)
    rop->scale(c1.scale());
  return true;
}

bool Evaluator::Impl::rotate_slots_inplace(Cipher*rop, RotKey const& rkey) const
{
  if (!rop)
    return false;
  if (rkey.galois == 0)
    return true;
  apply_galois(&(rop->bx), rkey.galois);
  apply_galois(&(rop->ax), rkey.galois);
  std::array<context::poly_t *, 2> result{&rop->bx, &rop->ax};
  std::array<context::poly_t *, 2>  input{&rop->bx, &rop->ax};
  apply_rotation_key(result, input, rkey);
  return true;
}

bool Evaluator::Impl::rotate_slots_inplace_non_ntt(Cipher *rop, RotKey const& rkey) const
{
  if (!rop)
    return false;
  if (rkey.galois == 0)
    return true;
  apply_galois_non_ntt(&(rop->bx), rkey.galois);
  apply_galois_non_ntt(&(rop->ax), rkey.galois);
  std::array<context::poly_t *, 2> result{&rop->bx, &rop->ax};
  std::array<context::poly_t *, 2>  input{&rop->bx, &rop->ax};
  apply_rotation_key_non_ntt(result, input, rkey);
  return true;
}

//! Input polynomials [d] should be in the ntt domain.
//! Resulting polynomials [rop] are also in the ntt domain.
void Evaluator::Impl::apply_rotation_key(
  std::array<context::poly_t *, 2> rop,
  std::array<context::poly_t *, 2> d,
  RotKey const &rotk) const 
{
  constexpr size_t L = context::nr_ctxt_moduli;
  constexpr size_t K = context::nr_sp_primes;
  if (!rop[0] || !rop[1] || !d[0] || !d[1])
    throw std::invalid_argument("Nullptr error.");
  const size_t Li = rop[0]->moduli_count();
  if (Li > L) 
    throw std::invalid_argument("Too many moduli.");
  assert(d[1]->moduli_count() == Li);
  context::poly_t d_rotk[2] = { context::poly_t(Li + K),
                                context::poly_t(Li + K) };

  //! convert back to power-basis for the mod up operation.
  d[1]->backward();
  bconv.approximated_mod_up(&d_rotk[0], (*d[1])); 
  //! forward to ntt domain for multiplication but now we have special moduli in d_rotk[0].
  for (size_t cm = 0; cm < Li + K; ++cm) {
    const size_t moduli_idx = cm < Li ? cm : context::index_sp_prime(cm - Li);
    yell::ntt<context::degree>::forward(d_rotk[0].ptr_at(cm), moduli_idx);
  }
  d_rotk[1] = d_rotk[0]; // copy

  //! multiply the d[1] to rotation key over the extened moduli.
  //! d[1] consists of Li (Li <= L) normal moduli and K special moduli,
  //! however the rotation key consists of L normal moduli and K special moduli.
  //! We need to carefully use the right moduli when doing the multiplication.
  for (size_t cm = 0; cm < Li + K; ++cm) {
    const size_t moduli_idx = cm < Li ? cm : context::index_sp_prime(cm - Li);
    auto  beta = rotk.beta.cptr_at(moduli_idx);
    auto alpha = rotk.alpha.cptr_at(moduli_idx);
    auto op0   = d_rotk[0].ptr_at(cm);
    auto op1   = d_rotk[1].ptr_at(cm);
    yell::ops::mulmod mulmod;
    for (size_t d = 0; d < context::degree; ++d) {
      mulmod.compute(*op0++, *beta++, moduli_idx);
      mulmod.compute(*op1++, *alpha++, moduli_idx);
    }
  }

  //! backward to power-basis for the mod down operation (with special moduli)
  for (size_t cm = 0; cm < Li + K; ++cm) {
    const size_t moduli_idx = cm < Li ? cm : context::index_sp_prime(cm - Li);
    yell::ntt<context::degree>::backward(d_rotk[0].ptr_at(cm), moduli_idx);
    yell::ntt<context::degree>::backward(d_rotk[1].ptr_at(cm), moduli_idx);
  }

  auto cpy_d0(*d[0]);
  //! mod down
  for (int i : {0, 1}) {
    bconv.approximated_mod_down(rop[i], d_rotk[i]); // rop[i] is now in the power-basis
    rop[i]->forward();
  }
  (*rop[0]) += cpy_d0;
}

//! Input polynomials [d] should be in the power domain.
//! Resulting polynomials [rop] are also in the power domain.
void Evaluator::Impl::apply_rotation_key_non_ntt(
  std::array<context::poly_t *, 2> rop,
  std::array<context::poly_t *, 2> d,
  RotKey const &rotk) const 
{
  constexpr size_t L = context::nr_ctxt_moduli;
  constexpr size_t K = context::nr_sp_primes;
  if (!rop[0] || !rop[1] || !d[0] || !d[1])
    throw std::invalid_argument("Nullptr error.");
  const size_t Li = rop[0]->moduli_count();
  if (Li > L) 
    throw std::invalid_argument("Too many moduli.");
  assert(d[1]->moduli_count() == Li);
  context::poly_t d_rotk[2] = { context::poly_t(Li + K),
                                context::poly_t(Li + K) };

  //! convert back to power-basis for the mod up operation.
  bconv.approximated_mod_up(&d_rotk[0], (*d[1])); 
  //! forward to ntt domain for multiplication but now we have special moduli in d_rotk[0].
  for (size_t cm = 0; cm < Li + K; ++cm) {
    const size_t moduli_idx = cm < Li ? cm : context::index_sp_prime(cm - Li);
    yell::ntt<context::degree>::forward(d_rotk[0].ptr_at(cm), moduli_idx);
  }
  d_rotk[1] = d_rotk[0]; // copy

  //! multiply the d[1] to rotation key over the extened moduli.
  //! d[1] consists of Li (Li <= L) normal moduli and K special moduli,
  //! however the rotation key consists of L normal moduli and K special moduli.
  //! We need to carefully use the right moduli when doing the multiplication.
  for (size_t cm = 0; cm < Li + K; ++cm) {
    const size_t moduli_idx = cm < Li ? cm : context::index_sp_prime(cm - Li);
    auto  beta = rotk.beta.cptr_at(moduli_idx);
    auto alpha = rotk.alpha.cptr_at(moduli_idx);
    auto op0   = d_rotk[0].ptr_at(cm);
    auto op1   = d_rotk[1].ptr_at(cm);
    yell::ops::mulmod mulmod;
    for (size_t d = 0; d < context::degree; ++d) {
      mulmod.compute(*op0++, *beta++, moduli_idx);
      mulmod.compute(*op1++, *alpha++, moduli_idx);
    }
  }

  //! backward to power-basis for the mod down operation (with special moduli)
  for (size_t cm = 0; cm < Li + K; ++cm) {
    const size_t moduli_idx = cm < Li ? cm : context::index_sp_prime(cm - Li);
    yell::ntt<context::degree>::backward(d_rotk[0].ptr_at(cm), moduli_idx);
    yell::ntt<context::degree>::backward(d_rotk[1].ptr_at(cm), moduli_idx);
  }

  auto cpy_d0(*d[0]); // d[0], rop[0] might point to the same object,
  //! mod down
  for (int i : {0, 1})
    bconv.approximated_mod_down(rop[i], d_rotk[i]); // rop[i] is now in the power-basis
  (*rop[0]) += cpy_d0;
}

bool Evaluator::Impl::multiply_inplace(
  Cipher *rop,
  Cipher const& c1) const
{
  if (!rop) return false;
  if (!can_mul_with(*rop, c1)) {
    std::cerr << "Can not multiply_inplace" << std::endl;
    return false;
  }
  using poly_t = context::poly_t;
  const size_t Li = rop->moduli_count();
  const double new_scale = rop->scale() / yell::params::P[Li - 1] * c1.scale();
  poly_t ct[3]{poly_t(Li), poly_t(Li), poly_t(Li)};

  Mul mul;
  mul.cross_mult(ct, *rop, c1);
  rop->bx = std::move(ct[0]);
  rop->ax = std::move(ct[1]);
  rop->mul_aux = std::make_shared<context::poly_t>(ct[2]);
  rop->scale(new_scale);
  return true;
}

bool Evaluator::Impl::relinerize_inplace(Cipher *rop, MultKey const& mkey) const
{
  if (!rop) return false;
  if (rop->canonical()) return false;
  using poly_t = context::poly_t;
  const size_t Li = rop->moduli_count();
  //! Relinear
  //! NOTE: rop here are in the NTT domain.
  Mul mul;
  std::array<poly_t *, 2> tuple{&rop->bx, &rop->ax};
  std::array<const poly_t *, 3> triple{&rop->bx, &rop->ax, rop->mul_aux.get()};
  mul.reLinearize(tuple, triple, mkey);
  //! Results from reLinearize_via_bit_decompose are in the NTT domain.
  Cipher tmp(Li - 1);
  std::array<poly_t *, 2> res{&tmp.bx, &tmp.ax};
  mul.rescale_ntt(res, tuple, bconv);
  rop->bx = std::move(tmp.bx);
  rop->ax = std::move(tmp.ax);
  rop->mul_aux = nullptr;
  return true;
}

void Evaluator::Impl::Mul::cross_mult(
  context::poly_t ans[3],
  Cipher const& c0,
  Cipher const& c1) const
{
  using T = yell::params::value_type;
  using gT = yell::params::gt_value_type;
  const size_t nmoduli = ans[0].moduli_count();
  for (size_t cm = 0; cm < nmoduli; ++cm) {
    auto   dst = ans[1].ptr_at(cm);
    auto b0_ptr = c0.bx.cptr_at(cm);
    auto a0_ptr = c0.ax.cptr_at(cm);
    auto b1_ptr = c1.bx.cptr_at(cm);
    auto a1_ptr = c1.ax.cptr_at(cm);

    for (size_t d = 0; d < context::degree; ++d) {
      gT v0 = (gT) (*b0_ptr++) * (*a1_ptr++);
      auto v1 = (gT) (*b1_ptr++) * (*a0_ptr++);
      v0 += v1;
      yell::ops::barret_reduction(&v0, cm);
      *dst++ = (T) v0;
    }
  }

  ans[0] = c0.bx; ans[0] *= c1.bx;
  ans[2] = c0.ax; ans[2] *= c1.ax;
}

void Evaluator::Impl::Mul::reLinearize(
  std::array<context::poly_t *, 2> &rop,
  std::array<const context::poly_t *, 3> ct,
  MultKey const& evk) const
{
  using T = yell::params::value_type;
  using gT = yell::params::gt_value_type;
  constexpr size_t degree= context::degree;
  constexpr size_t bytes = degree * sizeof(T);
  const size_t Li = ct[2]->moduli_count();
  auto ct2_power(*ct[2]);
  ct2_power.backward();
  //! evk.beta * ct[2], and evk.alpha * ct[2]
  std::array<T, degree> ct2;
  for (size_t i = 0; i < Li; ++i) {
    //! compute the i-th moduli of evk * ct[2]
    std::array<gT, degree> evk_ct2[2]{};
    for (size_t j = 0; j < Li; ++j) {
      std::memcpy(ct2.data(), ct2_power.cptr_at(j), bytes);
      yell::ntt<degree>::forward(ct2.data(), i);
      //! Lazy reduction for evk.beta * ct[2]
      //! Lazy reduction for evk.alpha * ct[2]
      //! Note: all these values are in the NTT domain.
      lazy_muladd(evk_ct2[0].begin(), 
                  evk.get_beta_at(j)->cptr_at(i),
                  ct2.cbegin(),
                  degree);
      lazy_muladd(evk_ct2[1].begin(), 
                  evk.get_alpha_at(j)->cptr_at(i),
                  ct2.cbegin(),
                  degree);
    }
    // Lazy reduction, ct[l] + evk.alpha * ct[l] for l \in {0, 1}
    for (int l : {0, 1}) {
      auto rop_ptr = rop[l]->ptr_at(i);
      auto ct_ptr = ct[l]->cptr_at(i);
      auto src_ptr = evk_ct2[l].cbegin();
      for (size_t d = 0; d < degree; ++d) {
        gT tmp = (*src_ptr++) + (*ct_ptr++);
        yell::ops::barret_reduction(&tmp, i);
        *rop_ptr++ = (T) tmp;
      }
    }
  }
}

void Evaluator::Impl::Mul::rescale_ntt(
  std::array<context::poly_t *, 2> rop,
  std::array<context::poly_t *, 2> ct,
  BaseConverter const& bconv) const
{
  for (int i : {0, 1}) {
    assert(rop[i] && ct[i]);
    assert(rop[i]->moduli_count() + 1 == ct[i]->moduli_count());
    ct[i]->backward();
    bconv.rescale(rop[i], *(ct[i]));
    rop[i]->forward();
  }
}

} // namespace fHE
