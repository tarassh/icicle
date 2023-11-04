//
//  ptx.metal
//  Test
//
//  Created by Taras Shchybovyk on 31.10.2023.
//

#include <metal_stdlib>
using namespace metal;

namespace ptx {
    // return x + y with uint32_t operands
    static inline uint32_t add(const uint32_t x, const uint32_t y) { return x + y; }

    // return x + y + carry with uint32_t operands
    static inline uint32_t addc(const uint32_t x, const uint32_t y, const uint32_t carry) { return x + y + carry; }

    // return x + y and carry out with uint32_t operands
    static inline uint32_t add_cc(const uint32_t x, const uint32_t y, thread uint32_t& carry)
    {
        uint32_t result;
        result = x + y;
        carry = x > result;
        return result;
    }

    // return x + y + carry and carry out  with uint32_t operands
    static inline uint32_t addc_cc(const uint32_t x, const uint32_t y, thread uint32_t& carry)
    {
        const uint32_t result = x + y + carry;
        carry = carry && x >= result || !carry && x > result;
        return result;
    }

    // return x - y with uint32_t operands
    static inline uint32_t sub(const uint32_t x, const uint32_t y) { return x - y; }

    //     return x - y - borrow with uint32_t operands
    static inline uint32_t subc(const uint32_t x, const uint32_t y, const uint32_t borrow) { return x - y - borrow; }

    //    return x - y and borrow out with uint32_t operands
    static inline uint32_t sub_cc(const uint32_t x, const uint32_t y, thread uint32_t& borrow)
    {
        uint32_t result;
        result = x - y;
        borrow = x < result;
        return result;
    }

    //    return x - y - borrow and borrow out with uint32_t operands
    static inline uint32_t subc_cc(const uint32_t x, const uint32_t y, thread uint32_t& borrow)
    {
        const uint32_t result = x - y - borrow;
        borrow = borrow && x <= result || !borrow && x < result;
        return result;
    }
    
    static inline uint32_t mul_lo(uint32_t x, uint32_t y) {
        return x * y;
    }
    
    static inline uint32_t mul_hi(uint32_t x, uint32_t y) {
        return mulhi(x, y);
    }
    
    static inline uint32_t mad_lo_cc(uint32_t x, uint32_t y, uint32_t z, thread uint32_t& carry) {
        uint64_t full_result = static_cast<uint64_t>(x) * static_cast<uint64_t>(y) + z;
        carry = (full_result > UINT32_MAX);
        return static_cast<uint32_t>(full_result);
    }
    
    static inline uint32_t mad_hi_cc(uint32_t x, uint32_t y, uint32_t z, thread uint32_t& carry) {
        uint64_t product = static_cast<uint64_t>(x) * static_cast<uint64_t>(y);
        uint32_t high_product = static_cast<uint>(product >> 32); // High 32 bits of the product
        uint64_t full_result = static_cast<uint64_t>(high_product) + z;
        carry = (full_result >> 32) != 0; // If there's anything in the high 32 bits, there was a carry out.
        return static_cast<uint32_t>(full_result);
    }
    
    static inline uint32_t madc_lo(uint32_t x, uint32_t y, uint32_t z, thread uint32_t &carry) {
        uint64_t full_product = static_cast<uint64_t>(x) * static_cast<uint64_t>(y);
        uint64_t full_result = (full_product & 0xFFFFFFFFull) + z + carry;
        carry = (full_product >> 32) != 0;
        return static_cast<uint32_t>(full_result);
    }
    
    static inline uint32_t madc_hi(uint32_t x, uint32_t y, uint32_t z, thread uint32_t &carry) {
        uint64_t full_product = static_cast<uint64_t>(x) * static_cast<uint64_t>(y);
        uint64_t high_product = full_product >> 32; // High 32 bits of the product
        uint64_t full_result = high_product + z + carry;
        carry = (full_result >> 32) != 0;
        return static_cast<uint>(full_result);
    }
    
    static inline uint32_t madc_lo_cc(uint32_t x, uint32_t y, uint32_t z, thread uint32_t &carry) {
        uint64_t full_product = static_cast<uint64_t>(x) * static_cast<uint64_t>(y);
        uint64_t low_product = full_product & 0xFFFFFFFFull; // Low 32 bits of the product
        uint64_t full_result = low_product + z + carry;
        carry = (full_result >> 32) != 0; // If the result is larger than 32 bits, we have a carry out
        return static_cast<uint32_t>(full_result);
    }
    
    static inline uint32_t madc_hi_cc(uint32_t x, uint32_t y, uint32_t z, thread uint32_t &carry) {
        uint64_t full_product = static_cast<uint64_t>(x) * static_cast<uint64_t>(y);
        uint32_t high_product = static_cast<uint>(full_product >> 32); // High 32 bits of the product
        uint64_t full_result = static_cast<uint64_t>(high_product) + z + carry;
        carry = (full_result >> 32) != 0; // Check if there's a carry-out
        return static_cast<uint32_t>(full_result);
    }
} // namespace ptx
