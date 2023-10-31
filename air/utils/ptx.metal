//
//  ptx.metal
//  Test
//
//  Created by Taras Shchybovyk on 31.10.2023.
//

#include <metal_stdlib>
using namespace metal;


namespace ptx {
    inline uint32_t sub_cc(const uint32_t x, const uint32_t y, bool carry_in, thread bool &carry_out) {
        uint32_t result = x - y - (carry_in ? 1 : 0);
        carry_out = (x < y) || ((x == y) && carry_in);
        return result;
    }
    
    inline uint64_t subc(const uint64_t x, const uint64_t y, thread bool &borrow) {
        uint64_t result = x - y;
        borrow = x < y;
        return result;
    }
    
    inline uint64_t subc_cc(const uint64_t x, const uint64_t y, thread bool &borrow) {
        uint64_t result = x - y;
        borrow = x < y; // Set borrow if x is less than y
        return result;
    }
    
    inline uint64_t add_cc(const uint64_t x, const uint64_t y, bool carry_in, thread bool& carry_out) {
        uint64_t result = x + y + (carry_in ? 1 : 0);
        carry_out = (result < x) || ((result == x) && carry_in);
        return result;
    }
    
    inline uint64_t addc(const uint64_t x, const uint64_t y, thread bool &carry) {
        uint64_t result = x + y + (carry ? 1 : 0);
        carry = (result < x) || (result == x && carry);  // Set carry if overflow, including previous carry
        return result;
    }
    
    inline uint64_t addc_cc(const uint64_t x, const uint64_t y, thread bool &carry) {
        uint64_t result = x + y + (carry ? 1 : 0);
        carry = result < x || (result == x && carry); // Set carry if overflow, including previous carry
        return result;
    }
} // namespace ptx
