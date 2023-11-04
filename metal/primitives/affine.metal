//
//  affine.metal
//  Test
//
//  Created by Taras Shchybovyk on 04.11.2023.
//

#pragma once
#include "field.metal"

template <class FF>
class Affine
{
public:
  FF x;
  FF y;

  static inline Affine neg(thread const Affine& point) { return {point.x, FF::neg(point.y)}; }

  friend inline bool operator==(thread const Affine& xs, thread const Affine& ys)
  {
    return (xs.x == ys.x) && (xs.y == ys.y);
  }
};

