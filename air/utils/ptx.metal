//
//  ptx.metal
//  Test
//
//  Created by Taras Shchybovyk on 31.10.2023.
//

#include <metal_stdlib>
using namespace metal;

namespace ptx {
  inline uint32_t sub_cc(const uint32_t x, const uint32_t y, bool carry_in, thread bool& carry_out)
  {
    uint32_t result = x - y - (carry_in ? 1 : 0);
    carry_out = (x < y) || ((x == y) && carry_in);
    return result;
  }

  inline uint64_t subc(const uint64_t x, const uint64_t y, thread bool& borrow)
  {
    uint64_t result = x - y;
    borrow = x < y;
    return result;
  }

  inline uint64_t subc_cc(const uint64_t x, const uint64_t y, thread bool& borrow)
  {
    uint64_t result = x - y;
    borrow = x < y; // Set borrow if x is less than y
    return result;
  }

  inline uint64_t add_cc(const uint64_t x, const uint64_t y, bool carry_in, thread bool& carry_out)
  {
    uint64_t result = x + y + (carry_in ? 1 : 0);
    carry_out = (result < x) || ((result == x) && carry_in);
    return result;
  }

  inline uint64_t addc(const uint64_t x, const uint64_t y, thread bool& carry)
  {
    uint64_t result = x + y + (carry ? 1 : 0);
    carry = (result < x) || (result == x && carry); // Set carry if overflow, including previous carry
    return result;
  }

  inline uint64_t addc_cc(const uint64_t x, const uint64_t y, thread bool& carry)
  {
    uint64_t result = x + y + (carry ? 1 : 0);
    carry = result < x || (result == x && carry); // Set carry if overflow, including previous carry
    return result;
  }

  inline uint64_t mul_lo(uint64_t x, uint64_t y)
  {
    uint32_t x_lo = x & 0xFFFFFFFF;
    uint32_t x_hi = x >> 32;
    uint32_t y_lo = y & 0xFFFFFFFF;
    uint32_t y_hi = y >> 32;

    uint64_t res_lo = uint64_t(x_lo) * y_lo;
    uint64_t res_mid_1 = uint64_t(x_lo) * y_hi;
    uint64_t res_mid_2 = uint64_t(x_hi) * y_lo;

    return res_lo + (res_mid_1 << 32) + (res_mid_2 << 32);
  }

  inline uint64_t mul_hi(uint64_t x, uint64_t y)
  {
    uint32_t x_lo = x & 0xFFFFFFFF;
    uint32_t x_hi = x >> 32;
    uint32_t y_lo = y & 0xFFFFFFFF;
    uint32_t y_hi = y >> 32;

    uint64_t hi = uint64_t(x_hi) * y_hi;
    uint64_t mid_1 = uint64_t(x_lo) * y_hi;
    uint64_t mid_2 = uint64_t(x_hi) * y_lo;

    // Carry bits from mid_1 and mid_2 can contribute to the high part
    uint64_t carry = (mid_1 >> 32) + (mid_2 >> 32);

    return hi + carry;
  }

  inline uint64_t mad_lo_cc(uint64_t x, uint64_t y, uint64_t z)
  {
    uint32_t x_lo = x & 0xFFFFFFFF;
    uint32_t x_hi = x >> 32;
    uint32_t y_lo = y & 0xFFFFFFFF;
    uint32_t y_hi = y >> 32;

    uint64_t product_lo = uint64_t(x_lo) * y_lo;
    uint64_t product_mid = uint64_t(x_lo) * y_hi + uint64_t(x_hi) * y_lo;
    uint64_t product_hi = uint64_t(x_hi) * y_hi;

    uint64_t sum = product_lo + (product_mid << 32) + z;

    // Overflow handling
    uint64_t carry = (product_lo > sum ? 1 : 0) + (product_mid >> 32) + (product_hi << 32);

    return sum;
  }

  inline uint64_t madc_hi_cc(uint64_t x, uint64_t y, uint64_t z)
  {
    uint32_t x_lo = x & 0xFFFFFFFF;
    uint32_t x_hi = x >> 32;
    uint32_t y_lo = y & 0xFFFFFFFF;
    uint32_t y_hi = y >> 32;

    uint64_t lo = uint64_t(x_lo) * y_lo;
    uint64_t mid_1 = uint64_t(x_lo) * y_hi;
    uint64_t mid_2 = uint64_t(x_hi) * y_lo;
    uint64_t hi = uint64_t(x_hi) * y_hi;

    // Calculate the carry from the lower bits
    uint64_t carry = (mid_1 >> 32) + (mid_2 >> 32);

    // Add the high parts and the carry
    uint64_t result = hi + (mid_1 >> 32) + (mid_2 >> 32) + (z >> 64);

    return result;
  }

  inline uint32_t madc_lo_cc(uint32_t x, uint32_t y, uint32_t z)
  {
    uint64_t product = static_cast<uint64_t>(x) * y;
    uint64_t sum = product + z;
    uint32_t result = static_cast<uint32_t>(sum & 0xFFFFFFFF);
    return result;
  }

  inline uint64_t madc_hi(uint64_t x, uint64_t y, uint64_t z)
  {
    uint32_t x_lo = x & 0xFFFFFFFF;
    uint32_t x_hi = x >> 32;
    uint32_t y_lo = y & 0xFFFFFFFF;
    uint32_t y_hi = y >> 32;

    uint64_t lo = uint64_t(x_lo) * y_lo;
    uint64_t mid_1 = uint64_t(x_lo) * y_hi;
    uint64_t mid_2 = uint64_t(x_hi) * y_lo;
    uint64_t hi = uint64_t(x_hi) * y_hi;

    // Calculate the carry from the lower bits
    uint64_t carry = (mid_1 >> 32) + (mid_2 >> 32);

    // Add the high parts, the carry and the z term
    uint64_t result = hi + carry + (z >> 64);

    return result;
  }
} // namespace ptx
