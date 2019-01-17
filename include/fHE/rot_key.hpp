#pragma once
#include "fHE/context.hpp"
#include <vector>
#include <memory>
#include <unordered_map>
namespace fHE {
struct SK; // forward declaration

struct RotKey {
private:
  static constexpr size_t degree = context::degree;
  static constexpr size_t L = context::nr_ctxt_moduli;
  static constexpr size_t K = context::nr_sp_primes;

  /*
     Multiply P (i.e., the product of all special primes).
     [x]_{q1, q2, ..., qL} --> [P * x]_{q1, q2, ..., qL}
   */
  void lift_normal_moduli(context::poly_t *op) const;

public:
  enum class Direction {
    LEFT,
    RIGHT
  };
  const size_t galois;
  context::poly_t beta; //! L + K moduli
  context::poly_t alpha;//! L + K moduli

  explicit RotKey(SK const &sk, uint16_t offset, Direction dir);

  RotKey(RotKey const&oth) = delete;

  RotKey(RotKey &&oth) = delete;

  RotKey& operator=(RotKey const&oth) = delete;

  static size_t rotation_offset_galois(int offset);
  ~RotKey();
};

struct RotKeySet {
private:
  static constexpr size_t logn = yell::static_log2<context::degree / 2>::value;
  using rotkey_ptr = std::shared_ptr<RotKey>;
  std::unordered_map<size_t, rotkey_ptr> rot_keys;

public:
  //! create log(N) (left) rotation key so that we can rotate with
  //! any offset in [0, N).
  explicit RotKeySet(SK const& sk);

  RotKeySet(RotKeySet const&oth) = delete;

  RotKeySet(RotKeySet &&oth) = delete;

  RotKeySet& operator=(RotKeySet const&oth) = delete;

  const RotKey* get(uint16_t offset) const;

  ~RotKeySet();
};

} // namespace fHE

