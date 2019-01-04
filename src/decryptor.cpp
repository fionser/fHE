#include "fHE/decryptor.hpp"
#include "fHE/cipher.hpp"
#include "fHE/keys.hpp"
#include "fHE/yell.hpp"
#include "fHE/encoder.hpp"
#include <algorithm>
namespace fHE {
Decryptor::Decryptor(std::shared_ptr<Encoder> encoder)
  : encoder(encoder) {}

bool Decryptor::decrypt(context::poly_t *rop, Cipher const& ctx, SK const& sk) const
{
  constexpr size_t degree = context::degree;
  using T = yell::params::value_type;
  using gT = yell::params::gt_value_type;
  const size_t nmoduli = ctx.moduli_count();
  if (nmoduli > sk.sx.moduli_count())
    throw std::runtime_error("Impossible.");
  if (!rop)
    return false;
  rop->resize_moduli_count(nmoduli);
  for (size_t cm = 0; cm < nmoduli; ++cm) {
    std::array<gT, degree> lazy;
    std::transform(ctx.bx.cptr_at(cm), ctx.bx.cptr_end(cm), lazy.data(),
                   [](T v) { return (gT) v; });
    //! m = bx + sx * ax
    lazy_muladd(lazy.data(), sk.sx.cptr_at(cm), ctx.ax.cptr_at(cm), degree);
    auto dst = rop->ptr_at(cm);
    for (auto &v : lazy) {
      yell::ops::barret_reduction(&v, cm);
      *dst++ = (T) v;
    }
  }
  rop->backward();
  return true;
}

bool Decryptor::decrypt(std::vector<std::complex<double>> *rop, 
                        Cipher const& ctx, 
                        SK const& sk) const
{
  if (!rop || !encoder)
    return false;
  context::poly_t plain(ctx.moduli_count());
  if (!decrypt(&plain, ctx, sk))
    return false;
  encoder->decode(rop, plain, ctx.scale());
  return true;
}

bool Decryptor::decrypt(std::vector<double> *rop, Cipher const& ctx, SK const& sk) const
{
  if (!rop || !encoder)
    return false;
  context::poly_t plain(ctx.moduli_count());
  if (!decrypt(&plain, ctx, sk))
    return false;
  encoder->decode(rop, plain, ctx.scale());
  return true;
}

} // namespace fHE
