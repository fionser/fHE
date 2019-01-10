#include "fHE/base_converter.hpp"
#include "fHE/yell.hpp"
#include <functional>
namespace fHE {
template <class T, size_t _1, size_t _2>
struct TupleArray {
  // a[_1][_2]
  using type = std::array<std::array<T, _2>, _1>;
};

template <class T, size_t _1, size_t _2, size_t _3>
struct TripleArray {
  // a[_1][_2][_3]
  using type = std::array<typename TupleArray<T, _2, _3>::type, _1>;
};

namespace internal {
template <typename T, size_t K, size_t L>
void do_approximated_basis_conversion(
  std::array<T, context::degree>  *rop,
  context::poly_t  const& op,
  std::array<T, K> const& hat_pi_invmod_pi,
  std::array<T, L> const& hat_pi_mod_qj,
  T const& P_mod_qj,
  std::function<size_t (size_t)> indexer,
  const size_t lcm);
} // namespace internal

struct BaseConverter::Impl {
public:
  using T = yell::params::value_type;
  /* number of special primes. */
  static constexpr size_t K = context::nr_sp_primes; 
  /* number of normal primes initially. */
  static constexpr size_t L = context::nr_ctxt_moduli; 
  static constexpr size_t degree = context::degree;

  /*
     Precomputed tables used for base conversions,
     mod_up and mod_down operations.
     We have L normal primes, q_1, q_2, ..., q_L.
     Also we have K special primes, p_1, p2, ..., p_K.
     Definition P := \prod_i p_i (the product of all special primes).
             ^p_i := P / p_i
                Q := \prod_j q_j (the product of all normal primes).
             ^q_j := Q / q_j
  */
  // Following three used to convert from special primes to a specific normal primes
  std::array<T, K>           hat_pi_invmod_pi; // (^p_i)^{-1} mod p_i
  TupleArray<T, L, K>::type  neg_hat_pi_mod_qj;// -(^p_i) mod q_j loop as [l][k]
  std::array<T, L>           P_invmod_qj;      // P^{-1} mod q_j
  std::array<T, L>           P_mod_qj;         // P mod q_j
  // Following three used to convert from normal primes to a specific special primes.
  // Since the number of normal primes in the chain would change, the first loop of
  // the above arrays indicates the number of normal primes.
  // Also, we do not convert from one normal prime, and thus we consider the case with at least 2 normal primes.
  TupleArray<T, L - 1, L>::type     hat_qj_invmod_qj; // (^q_j)^{-1} mod q_j
  TripleArray<T, L - 1, K, L>::type hat_qj_mod_pi;    // loop as [lcm][k][l]
  TupleArray<T, L - 1, K>::type     Q_mod_pi;         // Q mod p_i

  Impl();

  bool approximated_mod_down(context::poly_t *rop,
                             context::poly_t const& op) const;

  bool approximated_mod_up(context::poly_t *rop,
                           context::poly_t const& op) const;

  bool rescale(context::poly_t *rop, context::poly_t const& op) const;

private:
  void approx_convert_to_special_basis(std::array<T, degree> *rop,
                                       context::poly_t const& op,
                                       const size_t sp_moduli_index) const;

  void neg_approx_convert_to_normal_basis(std::array<T, degree> *rop,
                                          context::poly_t const& op,
                                          const size_t nrl_moduli_index) const;
};

BaseConverter::BaseConverter() : impl_(std::make_shared<Impl>()) { }

BaseConverter::~BaseConverter() {}

bool BaseConverter::approximated_mod_up(
  context::poly_t *rop,
  context::poly_t const& op) const
{
  if (!impl_) return false;
  return impl_->approximated_mod_up(rop, op);
}

bool BaseConverter::approximated_mod_down(context::poly_t *rop,
                           context::poly_t const& op) const
{
  if (!impl_) return false;
  return impl_->approximated_mod_down(rop, op);
}

bool BaseConverter::rescale(context::poly_t *rop, context::poly_t const& op) const
{
  if (!impl_) return false;
  return impl_->rescale(rop, op);
}

BaseConverter::Impl::Impl() {
  using namespace yell;
  ops::mulmod mulmod;
  std::vector<size_t> sp_moduli(K);
  std::iota(sp_moduli.begin(), sp_moduli.end(), context::index_sp_prime(0));
  // hat_pi_invmod_pi[k]
  for (size_t k = 0; k < K; ++k) { 
    auto moduli(sp_moduli);
    moduli.erase(moduli.begin() + k);
    size_t cm = context::index_sp_prime(k);
    auto hat_pi_mod_pi = product_of(moduli.begin(), moduli.end(), cm);
    hat_pi_invmod_pi.at(k) = math::inv_mod_prime(hat_pi_mod_pi, cm);
    // neg_hat_pi_mod_qj[l][k]
    for (size_t l = 0; l < L; ++l)
      neg_hat_pi_mod_qj[l][k] = yell::params::P[l] - product_of(moduli.begin(), moduli.end(), l);
  }

  // hat_qj_invmod_qj[l]
  std::memset(hat_qj_invmod_qj.data(), 0, sizeof(hat_qj_invmod_qj));
  for (size_t nr_nrml_primes = 2; nr_nrml_primes <= L; ++nr_nrml_primes) {
    std::vector<size_t> nrl_moduli(nr_nrml_primes);
    std::iota(nrl_moduli.begin(), nrl_moduli.end(), 0);
    size_t const l = nr_nrml_primes - 2;
    for (size_t j = 0; j < nr_nrml_primes; ++j) {
      std::vector<size_t> moduli(nrl_moduli);
      moduli.erase(moduli.begin() + j);
      auto hat_qj_qj = product_of(moduli.begin(), moduli.end(), j);
      hat_qj_invmod_qj.at(l).at(j) = math::inv_mod_prime(hat_qj_qj, j);
    }
  }

  // hat_qj_mod_pi[l][k][j]
  std::memset(hat_qj_mod_pi.data(), 0, sizeof(hat_qj_mod_pi));
  for (size_t nr_nrml_primes = 2; nr_nrml_primes <= L; ++nr_nrml_primes) {
    std::vector<size_t> nrl_moduli(nr_nrml_primes);
    std::iota(nrl_moduli.begin(), nrl_moduli.end(), 0);
    size_t const l = nr_nrml_primes - 2;
    for (size_t k = 0; k < K; ++k) {
      for (size_t j = 0; j < nr_nrml_primes; ++j) {
        std::vector<size_t> moduli{nrl_moduli};
        moduli.erase(moduli.begin() + j);
        const size_t kcm = context::index_sp_prime(k);
        hat_qj_mod_pi[l][k][j] = product_of(moduli.begin(), moduli.end(), kcm);
      }
    }
  }

  // P_invmod_qj
  for (size_t j = 0; j < L; ++j) {
    auto P = product_of(sp_moduli.begin(), sp_moduli.end(), j);
    P_mod_qj[j] = P;
    P_invmod_qj[j] = math::inv_mod_prime(P, j);
  }

  // Q_mod_pi
  std::memset(Q_mod_pi.data(), 0, sizeof(Q_mod_pi));
  for (size_t nr_nrml_primes = 2; nr_nrml_primes <= L; ++nr_nrml_primes) {
    size_t const lcm = nr_nrml_primes - 2;
    std::vector<size_t> nrl_moduli(nr_nrml_primes);
    std::iota(nrl_moduli.begin(), nrl_moduli.end(), 0);
    for (size_t k = 0; k < K; ++k) {
      size_t kcm = context::index_sp_prime(k);
      Q_mod_pi[lcm][k] = product_of(nrl_moduli.begin(), nrl_moduli.end(), kcm);
    }
  }
}

bool BaseConverter::Impl::approximated_mod_up(
  context::poly_t *rop,
  context::poly_t const& op) const
{
  if (!rop) return false;
  const size_t nmoduli = op.moduli_count();
  assert(nmoduli <= L && "Only support mod up from normal primes");
  for (size_t cm = 0; cm < nmoduli; ++cm)
    std::memcpy(rop->ptr_at(cm), op.cptr_at(cm), sizeof(T) * degree);
  for (size_t k = 0; k < K; ++k) {
    auto rop_ = yell::recast_as_array(*rop, nmoduli + k);
    approx_convert_to_special_basis(rop_, op, k);
  }
  return true;
}

bool BaseConverter::Impl::approximated_mod_down(
  context::poly_t *rop,
  context::poly_t const& op) const
{
  if (!rop) return false;
  assert(op.moduli_count() > K);
  const size_t L0 = rop->moduli_count();
  const size_t L1 = op.moduli_count();
  assert(L0 + K == L1);
  /* NOTE: point to the special moduli part. */
  context::poly_t sp_part(K);
  for (size_t k = 0; k < K; ++k)
    std::memcpy(sp_part.ptr_at(k), op.cptr_at(L0 + k), sizeof(T) * degree);

  // Convert basis from special moduli to L0 normal primes.
  for (size_t j = 0; j < L0; ++j) {
    auto rop_pointer = yell::recast_as_array(*rop, j);
    neg_approx_convert_to_normal_basis(rop_pointer, sp_part, j);
  }

  yell::ops::mulmod_shoup mulmod;
  for (size_t j = 0; j < L0; ++j) {
    auto dst = rop->ptr_at(j);
    auto src = op.cptr_at(j);
    auto end = op.cptr_end(j);
    auto shoup = yell::ops::shoupify(P_invmod_qj[j], j);
    while (src != end) {
      //! In math we compute (op - rop) * P_invmod_qj
      //! But the approx_convert_to_normal_basis compute the negative rop already
      //! Thus, we just need (lazy) addition here.
      auto temp = (*src++) + (*dst);
      *dst++ = mulmod(temp, P_invmod_qj[j], shoup, j);
    }
  }
  return true;
}

bool BaseConverter::Impl::rescale(context::poly_t *rop, context::poly_t const& op) const
{
  if (!rop) return false;
  using namespace yell;
  const size_t nmoduli = op.moduli_count();
  assert(nmoduli > 1 && "Already in the last moduli.");
  assert(rop->moduli_count() + 1 == nmoduli);
  yell::ops::mulmod_shoup mulmod;
  const auto last_prime = params::P[nmoduli - 1];
  // NOTE: for each moduli, sub the last moduli and then multiply qinv.
  for (size_t j = 0; j < nmoduli - 1; ++j) {
    //! take the last rns moduli with mod qj.
    auto last_mod = op.cptr_at(nmoduli - 1); 
    auto qinv = math::inv_mod_prime(last_prime, j);
    auto qinv_shoup = ops::shoupify(qinv, j);
    assert(qinv > 0);
    auto dst = rop->ptr_at(j);
    auto head = op.cptr_at(j);
    auto tail = last_mod;
    for (size_t d = 0; d < degree; ++d) {
      //! (head - tail) * qinv mod p_j, we save one reduction for subtraction.
      auto tmp = *head++ + params::P[j] - *tail++;
      *dst++ = mulmod(tmp, qinv, qinv_shoup, j);
    }
  }
  return true;
}

void BaseConverter::Impl::approx_convert_to_special_basis(
  std::array<T, degree> *rop,
  context::poly_t const& op,
  const size_t sp_moduli_index) const
{
  assert(sp_moduli_index < K);
  if (!rop) return;
  const size_t nmoduli = op.moduli_count();
  if (nmoduli <= 1)
    throw std::invalid_argument("nmoduli <= 1");
  const size_t lcm = nmoduli - 2; // since moduli > 1; and index starts from 0, so -2 here.
  const size_t kcm = context::index_sp_prime(sp_moduli_index);
  auto indexer = [](size_t l) { return l; };
  internal::do_approximated_basis_conversion(
    rop, op, 
    hat_qj_invmod_qj[lcm],
    hat_qj_mod_pi[lcm][sp_moduli_index],
    Q_mod_pi[lcm][sp_moduli_index],
    indexer, kcm);
}

void BaseConverter::Impl::neg_approx_convert_to_normal_basis(
  std::array<T, degree> *rop,
  context::poly_t const& op,
  const size_t nrl_moduli_index) const
{
  assert(nrl_moduli_index < L);
  if (!rop) return;
  auto indexer = [](size_t k) { return context::index_sp_prime(k); };
  internal::do_approximated_basis_conversion(
    rop, op, 
    hat_pi_invmod_pi,
    neg_hat_pi_mod_qj.at(nrl_moduli_index),
    P_mod_qj.at(nrl_moduli_index),
    indexer, nrl_moduli_index);
}

namespace internal {
template <typename T, size_t K, size_t L>
void do_approximated_basis_conversion(
  std::array<T, context::degree>  *rop,
  context::poly_t  const& op,
  std::array<T, K> const& hat_pi_invmod_pi,
  std::array<T, L> const& hat_pi_mod_qj,
  T const& P_mod_qj,
  std::function<size_t (size_t)> indexer,
  const size_t lcm) 
{
  if (!rop) return;
  const size_t nmoduli = op.moduli_count();
  assert(nmoduli <= K);

  using gT = yell::params::gt_value_type;
  yell::ops::mulmod_shoup mulmod;
  std::array<gT, context::degree> y{};
  for (size_t k = 0; k < nmoduli; ++k) {
    const size_t kcm = indexer(k);
    const T hpi_inv_pi = hat_pi_invmod_pi[k];
    const T shoup = yell::ops::shoupify(hpi_inv_pi, kcm);
    const T hpi_qj = hat_pi_mod_qj[k];
    auto x = op.cptr_at(k);
    for (auto &yd : y) {
      //! x * hat_qj_invmod_qj mod qj
      T t0 = mulmod(*x++, hpi_inv_pi, shoup, kcm);
      //! lazy reduction
      yd += (gT) t0 * hpi_qj;
    }
  }

  for (size_t d = 0; d < context::degree; ++d) {
    yell::ops::barret_reduction(&y[d], lcm);
    (*rop)[d] = (T) y[d];
  }
}
} // namespace internal

} // namespace fHE
