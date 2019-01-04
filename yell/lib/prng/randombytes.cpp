#include "yell/prng/randombytes.h"
#if defined(__APPLE__) || defined(__linux__)
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#else
#include "../randutils.hpp"
#endif
namespace yell {
#if defined(__APPLE__) || defined(__linux__)
static int fd = -1;

void randombytes(unsigned char *x, unsigned long long xlen) {
  if (fd == -1) {
    for (;;) {
      fd = open("/dev/urandom", O_RDONLY);
      if (fd != -1) break;
      sleep(1);
    }
  }

  while (xlen > 0) {
    int i = (xlen < 1048576) ? xlen : 1048576;
    i = read(fd, x, i);

    if (i < 1) {
      sleep(1);
      continue;
    }

    x += i;
    xlen -= i;
  }
}
#else
#warning "Not using strong random source."
void randombytes(unsigned char *x, unsigned long long xlen) {
  randutils::mt19937_rng rng;
  rng.generate(x, x + xlen);
}
#endif
} // namespace nfl
