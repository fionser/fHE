#pragma once
#include "fHE/context.hpp"
#include <memory>
namespace fHE {
struct BaseConverter {
private:
  struct Impl;
  std::shared_ptr<Impl> impl_;

public:
  using T = yell::params::value_type;
  static constexpr size_t degree = context::degree;
  BaseConverter();

  ~BaseConverter();

  //! adding special primes
  bool approximated_mod_up(context::poly_t *rop,
                           context::poly_t const& op) const;

  bool exact_mod_up(context::poly_t *rop, context::poly_t const& op) const;

  // The special primes part is in the power-basis
  // The normal primes part is in the NTT-basis
  bool approximated_mod_down(context::poly_t *rop,
                             context::poly_t const& op) const;

  // The last moduli of [op] is in the power-basis, and other moduli are in the NTT-basis
  bool rescale(context::poly_t *rop, context::poly_t const& op) const;

  void approx_convert_to_special_basis(std::array<T, degree> *rop,
                                       context::poly_t const& op,
                                       const size_t sp_moduli_index) const;

  void convert_to_special_basis(std::array<T, degree> *rop,
                                context::poly_t const& op,
                                const size_t sp_moduli_index) const;

  void neg_approx_convert_to_normal_basis(std::array<T, degree> *rop,
                                          context::poly_t const& op,
                                          const size_t nrl_moduli_index) const;

};
} // namespace fHE
