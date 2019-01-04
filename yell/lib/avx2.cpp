#ifdef YELL_USE_AVX2
#include "yell/avx2.h"
#include <immintrin.h>

__m256i avx_mm256_mul64_hi(__m256i const& A, __m256i const& B)
{
  //! A := A_hi * 2^32 + A_lo;
  //! B := B_hi * 2^32 + B_lo; 
  //! Compute (A * B) >> 64
  static const __m256i avx_0F_mask = _mm256_set1_epi64x(0x00000000FFFFFFFFULL);
  //! exchange the high and low 32-bits
  __m256i Ahi = _mm256_shuffle_epi32(A, 1 | (0 << 2) | (3 << 4) | (2 << 6));
  __m256i Bhi = _mm256_shuffle_epi32(B, 1 | (0 << 2) | (3 << 4) | (2 << 6));

  __m256i AhiBhi = _mm256_mul_epu32(Ahi, Bhi);
  __m256i AhiBlo = _mm256_mul_epu32(Ahi, B);
  __m256i AloBhi = _mm256_mul_epu32(A, Bhi);

  __m256i AloBlo_hi = _mm256_srli_epi64(_mm256_mul_epu32(A,  B), 32);
  __m256i AhiBlo_lo = _mm256_and_si256(avx_0F_mask, AhiBlo);
  __m256i AloBhi_lo = _mm256_and_si256(avx_0F_mask, AloBhi);
  __m256i carry = _mm256_srli_epi64(_mm256_add_epi64(AloBlo_hi, _mm256_add_epi64(AhiBlo_lo, AloBhi_lo)), 32);

  __m256i AhiBlo_hi = _mm256_srli_epi64(AhiBlo, 32);
  __m256i AloBhi_hi = _mm256_srli_epi64(AloBhi, 32);

  return _mm256_add_epi64(carry, _mm256_add_epi64(AhiBhi, _mm256_add_epi64(AhiBlo_hi, AloBhi_hi)));
}

__m256i avx_mm256_mul64_lo(__m256i const& A, __m256i const& B)
{
  //! A := A_hi * 2^32 + A_lo;
  //! B := B_hi * 2^32 + B_lo; 
  //! Compute (A * B) bmod 2^64
  __m256i Ahi = _mm256_shuffle_epi32(A, 1 | (0 << 2) | (3 << 4) | (2 << 6));
  __m256i Bhi = _mm256_shuffle_epi32(B, 1 | (0 << 2) | (3 << 4) | (2 << 6));
  __m256i AloBlo = _mm256_mul_epu32(A, B);
  __m256i AhiBlo = _mm256_slli_epi64(_mm256_mul_epu32(Ahi, B), 32);
  __m256i AloBhi = _mm256_slli_epi64(_mm256_mul_epu32(A, Bhi), 32);
  return _mm256_add_epi64(AloBlo, _mm256_add_epi64(AhiBlo, AloBhi));
}

__m256i avx_mm256_mul64_hi(__m256i const& A, __m256i const& Ahi, __m256i const& B)
{
  //! A := A_hi * 2^32 + A_lo;
  //! B := B_hi * 2^32 + B_lo; 
  //! Compute (A * B) >> 64
  static const __m256i avx_0F_mask = _mm256_set1_epi64x(0x00000000FFFFFFFFULL);
  //! exchange the high and low 32-bits
  __m256i Bhi = _mm256_shuffle_epi32(B, 1 | (0 << 2) | (3 << 4) | (2 << 6));

  __m256i AhiBhi = _mm256_mul_epu32(Ahi, Bhi);
  __m256i AhiBlo = _mm256_mul_epu32(Ahi, B);
  __m256i AloBhi = _mm256_mul_epu32(A, Bhi);

  __m256i AloBlo_hi = _mm256_srli_epi64(_mm256_mul_epu32(A,  B), 32);
  __m256i AhiBlo_lo = _mm256_and_si256(avx_0F_mask, AhiBlo);
  __m256i AloBhi_lo = _mm256_and_si256(avx_0F_mask, AloBhi);
  __m256i carry = _mm256_srli_epi64(_mm256_add_epi64(AloBlo_hi, _mm256_add_epi64(AhiBlo_lo, AloBhi_lo)), 32);

  __m256i AhiBlo_hi = _mm256_srli_epi64(AhiBlo, 32);
  __m256i AloBhi_hi = _mm256_srli_epi64(AloBhi, 32);

  return _mm256_add_epi64(carry, _mm256_add_epi64(AhiBhi, _mm256_add_epi64(AhiBlo_hi, AloBhi_hi)));
}

__m256i avx_mm256_mul64_lo(__m256i const& A, __m256i const& Ahi, __m256i const& B)
{
  //! A := A_hi * 2^32 + A_lo;
  //! B := B_hi * 2^32 + B_lo; 
  //! Compute (A * B) bmod 2^64
  __m256i Bhi = _mm256_shuffle_epi32(B, 1 | (0 << 2) | (3 << 4) | (2 << 6));
  __m256i AloBlo = _mm256_mul_epu32(A, B);
  __m256i AhiBlo = _mm256_slli_epi64(_mm256_mul_epu32(Ahi, B), 32);
  __m256i AloBhi = _mm256_slli_epi64(_mm256_mul_epu32(A, Bhi), 32);
  return _mm256_add_epi64(AloBlo, _mm256_add_epi64(AhiBlo, AloBhi));
}
#endif // YELL_USE_AVX2
