#pragma once
#include "fHE/context.hpp"
#include <memory>
namespace fHE {
struct Cipher;
struct MultKey;
struct RotKey;
struct MixedRotKey;
template<class T> struct RotKeySet;
struct Encoder;

struct Evaluator {
private:
  struct Impl;
  std::shared_ptr<Impl> impl_;// private implementation
  std::shared_ptr<Encoder> encoder;

public:
  explicit Evaluator(std::shared_ptr<Encoder> encoder);
  ~Evaluator();

  bool add(Cipher *rop, Cipher const& c1) const;
  bool multiply(Cipher *rop, Cipher const& c1) const;
  bool multiply_plain(Cipher *rop, context::poly_t const& p, const double scale = -1) const;
  bool relinearize(Cipher *rop, MultKey const& mkey) const;
  bool multiply_and_relin(Cipher *rop, Cipher const& c1, MultKey const& mkey) const;
  //! The rotation offset is specified in the RotKey.
  bool rotate_slots(Cipher *rop, RotKey const& rkey) const;
  bool rotate_slots_non_ntt(Cipher *rop, RotKey const& rkey) const;
  bool rotate_slots(Cipher *rop, MixedRotKey const& rkey) const;
  //! Rotation with arbitary offset
  bool rotate_slots(Cipher *rop, int offset, RotKeySet<RotKey> const& rkeys) const;
  bool rotate_slots(Cipher *rop, int offset, RotKeySet<MixedRotKey> const& rkeys) const;
  template <class RotKey>
  bool rotate_slots_non_ntt(Cipher *rop, int offset, RotKeySet<RotKey> const& rkeys) const;
  //! Replicate the specific position of slots, using logN rotations and 
  //! a single masking.
  template <class RotKey>
  bool replicate(Cipher *rop, size_t pos, RotKeySet<RotKey> const &rkeys) const;
};
} // namespace fHE
