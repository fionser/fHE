#pragma once
#include "fHE/context.hpp"
#include <memory>
namespace fHE {
struct BaseConverter {
private:
  struct Impl;
  std::shared_ptr<Impl> impl_;
public:
  BaseConverter();

  ~BaseConverter();

  bool approximated_mod_up(context::poly_t *rop,
                           context::poly_t const& op) const;

  bool approximated_mod_down(context::poly_t *rop,
                             context::poly_t const& op) const;

  bool rescale(context::poly_t *rop, context::poly_t const& op) const;
};
} // namespace fHE
