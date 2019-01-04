#include "fHE/encryptor.hpp"
#include "fHE/cipher.hpp"
#include "fHE/keys.hpp"
#include "fHE/encoder.hpp"
namespace fHE {
Encryptor::Encryptor(std::shared_ptr<Encoder> encoder)
  : encoder(encoder) {}

bool Encryptor::encrypt(Cipher *rop, context::poly_t const& msg, PK const& pk) const
{
  if (!rop)
    return false;
  const size_t msg_moduli = msg.moduli_count();
  const size_t ctx_moduli = rop->moduli_count();
  context::poly_t u(ctx_moduli, yell::ZO_dist());
  u.forward();
  rop->bx.set(context::gauss_struct(&context::fg_prng)); // e0
  //! Add msg to e0, since addition is also fine in the power-basis, we save one NTT here.
  if (msg_moduli == ctx_moduli) {
    (rop->bx) += msg; // e0 + msg
  } else {
    yell::ops::addmod addmod;
    for (size_t cm = 0; cm < ctx_moduli; ++cm) {
      auto msg_ptr = msg.cptr_at(0); // use the first moduli only
      auto bx_ptr = rop->bx.ptr_at(cm);
      for (size_t d = 0; d < context::degree; ++d)
        addmod.compute(*bx_ptr++, *msg_ptr++, cm);
    }
  }
  rop->bx.forward_lazy(); // lazy reduction
  rop->bx.add_product_of(u, pk.bx); // ctx.bx = u * pk.bx + e0 + msg

  rop->ax.set(context::gauss_struct(&context::fg_prng)); // e1
  rop->ax.forward_lazy(); // lazy reduction
  rop->ax.add_product_of(u, pk.ax); // ctx.ax = u * pk.ax + e1
  return true;
}

bool Encryptor::encrypt(Cipher *rop, std::vector<double> const& values, PK const& pk) const
{
  if (!rop || !encoder)
    return false;
  context::poly_t plain(context::nr_ctxt_moduli);
  encoder->encode(&plain, values, context::encoder_scale);
  encrypt(rop, plain, pk);
  rop->scale(context::encoder_scale);
  return true;
}

bool Encryptor::encrypt(Cipher *rop, 
                        std::vector<std::complex<double>> const& values, 
                        PK const& pk) const
{
  if (!rop || !encoder)
    return false;
  context::poly_t plain(context::nr_ctxt_moduli);
  encoder->encode(&plain, values, context::encoder_scale);
  encrypt(rop, plain, pk);
  rop->scale(context::encoder_scale);
  return true;
}

} // namespace fHE
