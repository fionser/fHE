#pragma once
#include "fHE/context.hpp"
#include <vector>
#include <complex>
#include <memory>
namespace fHE {
struct SK;
struct Cipher;
struct Encoder;

struct Decryptor {
  explicit Decryptor(std::shared_ptr<Encoder> encoder = nullptr);

  bool decrypt(context::poly_t *rop, Cipher const& ctx, SK const& sk) const;

  bool decrypt(std::vector<double> *rop, Cipher const& ctx, SK const& sk) const;
  bool decrypt(std::vector<std::complex<double>> *rop, Cipher const& ctx, SK const& sk) const;

  std::shared_ptr<Encoder> encoder;
};
} // namespace fHE
