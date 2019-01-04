#pragma once
#include "fHE/context.hpp"
#include <memory>
namespace fHE {
struct Cipher;
struct MultKey;
struct RotKey;

struct Evaluator {
private:
  struct Impl;
  std::shared_ptr<Impl> impl_;// private implementation
public:
  Evaluator();
  ~Evaluator();
  bool add_inplace(Cipher *rop, Cipher const& c1) const;
  bool multiply_inplace(Cipher *rop, Cipher const& c1, MultKey const& mkey) const;
  bool rotate_slots_inplace(Cipher*rop, RotKey const& rkey) const;
};
} // namespace fHE
