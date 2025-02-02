//
//  projective.metal
//  Test
//
//  Created by Taras Shchybovyk on 04.11.2023.
//

#pragma once

#include "affine.metal"

template <typename FF, class SCALAR_FF, thread const FF& B_VALUE, thread const FF& GENERATOR_X, thread const FF& GENERATOR_Y>
class Projective
{
  friend Affine<FF>;

public:
  FF x;
  FF y;
  FF z;

  static inline Projective zero() { return {FF::zero(), FF::one(), FF::zero()}; }

  static inline Affine<FF> to_affine(thread const Projective& point)
  {
    FF denom = FF::inverse(point.z);
    return {point.x * denom, point.y * denom};
  }

  static inline Projective from_affine(device const Affine<FF>& point) { return {point.x, point.y, FF::one()}; }

  static inline Projective generator() { return {GENERATOR_X, GENERATOR_Y, FF::one()}; }

  static inline Projective neg(thread const Projective& point) { return {point.x, FF::neg(point.y), point.z}; }

  friend inline Projective operator+(Projective p1, device const Projective& p2)
  {
    const FF X1 = p1.x;                                                                //                   < 2
    const FF Y1 = p1.y;                                                                //                   < 2
    const FF Z1 = p1.z;                                                                //                   < 2
    const FF X2 = p2.x;                                                                //                   < 2
    const FF Y2 = p2.y;                                                                //                   < 2
    const FF Z2 = p2.z;                                                                //                   < 2
    const FF t00 = X1 * X2;                                                            // t00 ← X1 · X2     < 2
    const FF t01 = Y1 * Y2;                                                            // t01 ← Y1 · Y2     < 2
    const FF t02 = Z1 * Z2;                                                            // t02 ← Z1 · Z2     < 2
    const FF t03 = X1 + Y1;                                                            // t03 ← X1 + Y1     < 4
    const FF t04 = X2 + Y2;                                                            // t04 ← X2 + Y2     < 4
    const FF t05 = t03 * t04;                                                          // t03 ← t03 · t04   < 3
    const FF t06 = t00 + t01;                                                          // t06 ← t00 + t01   < 4
    const FF t07 = t05 - t06;                                                          // t05 ← t05 − t06   < 2
    const FF t08 = Y1 + Z1;                                                            // t08 ← Y1 + Z1     < 4
    const FF t09 = Y2 + Z2;                                                            // t09 ← Y2 + Z2     < 4
    const FF t10 = t08 * t09;                                                          // t10 ← t08 · t09   < 3
    const FF t11 = t01 + t02;                                                          // t11 ← t01 + t02   < 4
    const FF t12 = t10 - t11;                                                          // t12 ← t10 − t11   < 2
    const FF t13 = X1 + Z1;                                                            // t13 ← X1 + Z1     < 4
    const FF t14 = X2 + Z2;                                                            // t14 ← X2 + Z2     < 4
    const FF t15 = t13 * t14;                                                          // t15 ← t13 · t14   < 3
    const FF t16 = t00 + t02;                                                          // t16 ← t00 + t02   < 4
    const FF t17 = t15 - t16;                                                          // t17 ← t15 − t16   < 2
    const FF t18 = t00 + t00;                                                          // t18 ← t00 + t00   < 2
    const FF t19 = t18 + t00;                                                          // t19 ← t18 + t00   < 2
    const FF t20 = FF::template mul_unsigned<3>(FF::template mul_const<B_VALUE>(t02)); // t20 ← b3 · t02    < 2
    const FF t21 = t01 + t20;                                                          // t21 ← t01 + t20   < 2
    const FF t22 = t01 - t20;                                                          // t22 ← t01 − t20   < 2
    const FF t23 = FF::template mul_unsigned<3>(FF::template mul_const<B_VALUE>(t17)); // t23 ← b3 · t17    < 2
    const auto t24 = FF::mul_wide(t12, t23);                                           // t24 ← t12 · t23   < 2
    const auto t25 = FF::mul_wide(t07, t22);                                           // t25 ← t07 · t22   < 2
    const FF X3 = FF::reduce(t25 - t24);                                               // X3 ← t25 − t24    < 2
    const auto t27 = FF::mul_wide(t23, t19);                                           // t27 ← t23 · t19   < 2
    const auto t28 = FF::mul_wide(t22, t21);                                           // t28 ← t22 · t21   < 2
    const FF Y3 = FF::reduce(t28 + t27);                                               // Y3 ← t28 + t27    < 2
    const auto t30 = FF::mul_wide(t19, t07);                                           // t30 ← t19 · t07   < 2
    const auto t31 = FF::mul_wide(t21, t12);                                           // t31 ← t21 · t12   < 2
    const FF Z3 = FF::reduce(t31 + t30);                                               // Z3 ← t31 + t30    < 2
    return {X3, Y3, Z3};
  }

  friend inline Projective operator-(Projective p1, thread const Projective& p2) { return p1 + neg(p2); }

  friend inline Projective operator+(Projective p1, thread const Affine<FF>& p2)
  {
    const FF X1 = p1.x;                                                                //                   < 2
    const FF Y1 = p1.y;                                                                //                   < 2
    const FF Z1 = p1.z;                                                                //                   < 2
    const FF X2 = p2.x;                                                                //                   < 2
    const FF Y2 = p2.y;                                                                //                   < 2
    const FF t00 = X1 * X2;                                                            // t00 ← X1 · X2     < 2
    const FF t01 = Y1 * Y2;                                                            // t01 ← Y1 · Y2     < 2
    const FF t02 = Z1;                                                                 // t02 ← Z1          < 2
    const FF t03 = X1 + Y1;                                                            // t03 ← X1 + Y1     < 4
    const FF t04 = X2 + Y2;                                                            // t04 ← X2 + Y2     < 4
    const FF t05 = t03 * t04;                                                          // t03 ← t03 · t04   < 3
    const FF t06 = t00 + t01;                                                          // t06 ← t00 + t01   < 4
    const FF t07 = t05 - t06;                                                          // t05 ← t05 − t06   < 2
    const FF t08 = Y1 + Z1;                                                            // t08 ← Y1 + Z1     < 4
    const FF t09 = Y2 + FF::one();                                                     // t09 ← Y2 + 1      < 4
    const FF t10 = t08 * t09;                                                          // t10 ← t08 · t09   < 3
    const FF t11 = t01 + t02;                                                          // t11 ← t01 + t02   < 4
    const FF t12 = t10 - t11;                                                          // t12 ← t10 − t11   < 2
    const FF t13 = X1 + Z1;                                                            // t13 ← X1 + Z1     < 4
    const FF t14 = X2 + FF::one();                                                     // t14 ← X2 + 1      < 4
    const FF t15 = t13 * t14;                                                          // t15 ← t13 · t14   < 3
    const FF t16 = t00 + t02;                                                          // t16 ← t00 + t02   < 4
    const FF t17 = t15 - t16;                                                          // t17 ← t15 − t16   < 2
    const FF t18 = t00 + t00;                                                          // t18 ← t00 + t00   < 2
    const FF t19 = t18 + t00;                                                          // t19 ← t18 + t00   < 2
    const FF t20 = FF::template mul_unsigned<3>(FF::template mul_const<B_VALUE>(t02)); // t20 ← b3 · t02    < 2
    const FF t21 = t01 + t20;                                                          // t21 ← t01 + t20   < 2
    const FF t22 = t01 - t20;                                                          // t22 ← t01 − t20   < 2
    const FF t23 = FF::template mul_unsigned<3>(FF::template mul_const<B_VALUE>(t17)); // t23 ← b3 · t17    < 2
    const auto t24 = FF::mul_wide(t12, t23);                                           // t24 ← t12 · t23   < 2
    const auto t25 = FF::mul_wide(t07, t22);                                           // t25 ← t07 · t22   < 2
    const FF X3 = FF::reduce(t25 - t24);                                               // X3 ← t25 − t24    < 2
    const auto t27 = FF::mul_wide(t23, t19);                                           // t27 ← t23 · t19   < 2
    const auto t28 = FF::mul_wide(t22, t21);                                           // t28 ← t22 · t21   < 2
    const FF Y3 = FF::reduce(t28 + t27);                                               // Y3 ← t28 + t27    < 2
    const auto t30 = FF::mul_wide(t19, t07);                                           // t30 ← t19 · t07   < 2
    const auto t31 = FF::mul_wide(t21, t12);                                           // t31 ← t21 · t12   < 2
    const FF Z3 = FF::reduce(t31 + t30);                                               // Z3 ← t31 + t30    < 2
    return {X3, Y3, Z3};
  }

  friend inline Projective operator-(Projective p1, thread const Affine<FF>& p2)
  {
    return p1 + Affine<FF>::neg(p2);
  }

  friend inline Projective operator*(SCALAR_FF scalar, device const Projective& point)
  {
    Projective res = zero();
#pragma unroll
    for (uint i = 0; i < SCALAR_FF::NBITS; i++) {
      if (i > 0) { res = res + res; }
      if (scalar.get_scalar_digit(SCALAR_FF::NBITS - i - 1, 1)) { res = res + point; }
    }
    return res;
  }

  friend inline bool operator==(thread const Projective& p1, thread const Projective& p2)
  {
    return (p1.x * p2.z == p2.x * p1.z) && (p1.y * p2.z == p2.y * p1.z);
  }

  static inline bool is_zero(thread const Projective& point)
  {
    return point.x == FF::zero() && point.y != FF::zero() && point.z == FF::zero();
  }

  static inline bool is_on_curve(thread const Projective& point)
  {
    if (is_zero(point)) return true;
    bool eq_holds =
      (FF::template mul_const<B_VALUE>(FF::sqr(point.z) * point.z) + FF::sqr(point.x) * point.x ==
       point.z * FF::sqr(point.y));
    return point.z != FF::zero() && eq_holds;
  }

  static inline Projective rand_host()
  {
    SCALAR_FF rand_scalar = SCALAR_FF::rand_host();
    return rand_scalar * generator();
  }
};
