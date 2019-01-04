#include "yell/poly.hpp"
constexpr size_t Deg = 2048;
constexpr size_t NdM = 8;

namespace global {
  using gauss_struct = yell::gaussian<uint16_t, 2>;
  using gauss_t = yell::FastGaussianNoise<uint16_t, 2>;
  gauss_t fg_prng(3.2, 128, 1 << 14);
}

struct SK {

  explicit SK() : sx(NdM, yell::hwt_dist(64)) {
    sx.forward();
  }

  yell::poly<Deg> sx;
};

struct PK {

  explicit PK(SK &sk) : bx(NdM), ax(NdM, yell::uniform{}) {
    bx.set(global::gauss_struct(&global::fg_prng));
    bx.forward();
    bx.add_product_of(ax, sk.sx);
    bx.negate();
  }

  yell::poly<Deg> bx, ax;
};

struct CTXT {
  yell::poly<Deg> bx, ax;

  explicit CTXT(yell::poly<Deg> const& msg, PK const& pk)
    :bx(NdM), ax(NdM)
  {
    yell::poly<Deg> u(NdM, yell::ZO_dist{});
    u.forward();
    bx.set(global::gauss_struct(&global::fg_prng));
    bx += msg;
    bx.forward_lazy();
    bx.add_product_of(u, pk.bx);

    ax.set(global::gauss_struct(&global::fg_prng));
    ax.forward_lazy();
    ax.add_product_of(u, pk.ax);
  }
};

void decrypt(yell::poly<Deg> *rop, CTXT const& ctx, SK const& sk)
{
  if (!rop) return;
  (*rop) = ctx.bx;
  rop->add_product_of(sk.sx, ctx.ax);
  rop->backward();
}

int main() {
  SK sk;
  PK pk(sk);
  yell::poly<Deg> msg(NdM, yell::uniform{});

  CTXT ctx(msg, pk);
  yell::poly<Deg> plain(NdM);
  auto coeffs = plain.poly2mpz();

  decrypt(&plain, ctx, sk);

  std::cout << msg(0, 0) << " " << plain(0, 0) << "\n";
  return 0;
}
