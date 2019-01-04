#include <cassert>
namespace yell { namespace math {

uint32_t reverse_bits(uint32_t operand)
{
  operand = (((operand & 0xaaaaaaaa) >> 1) | ((operand & 0x55555555) << 1));
  operand = (((operand & 0xcccccccc) >> 2) | ((operand & 0x33333333) << 2));
  operand = (((operand & 0xf0f0f0f0) >> 4) | ((operand & 0x0f0f0f0f) << 4));
  operand = (((operand & 0xff00ff00) >> 8) | ((operand & 0x00ff00ff) << 8));
  return((operand >> 16) | (operand << 16));
}

/* The operand is less than 32 */

uint32_t reverse_bits(uint32_t operand, int32_t bit_count)
{
  assert(bit_count < 32);
  return (uint32_t) (((uint64_t) reverse_bits(operand)) >> (32 - bit_count));
}
}} // namespace yell::math
