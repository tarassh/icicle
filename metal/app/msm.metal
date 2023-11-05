//
//  msm.metal
//  Test
//
//  Created by Taras Shchybovyk on 31.10.2023.
//
//
#pragma once

#include <metal_stdlib>
using namespace metal;


template <typename P, typename S>
void ssm_kernel(device S* scalars,
                       device P* points,
                       device P* results,
                       const uint N,
                       uint tid) {
    if (tid < N) {
        results[tid] = scalars[tid] * points[tid];
    }
}

template <typename P, typename A>
void to_proj_kernel(device A* affine_points,
                    device P* proj_points,
                    const uint N,
                    uint tid) {
    if (tid < N) {
        proj_points[tid] = P::from_affine(affine_points[tid]);
    }
}

template <typename P>
void sum_reduction_kernel(device P* v,
                          device P* v_r,
                          device P* shared,
                          uint group_size,
                          uint tid,
                          uint tidx) {
    // Perform first level of reduction:
    // Each thread loads one element from global to shared mem
    
    // Load element and perform first level of reduction:
    // Each thread loads one element from global to shared mem
    shared[tidx] = v[tid];
    threadgroup_barrier(mem_flags::mem_threadgroup);
    
    // Start at 1/2 block stride and divide by two each iteration
    for (unsigned s = group_size / 2; s > 0; s >>= 1) {
        // Each thread does work unless it is further than the stride
        if (tidx < s) {
            shared[tidx] = shared[tidx] + shared[tidx + s];
        }
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
    
    // Let the thread 0 for this block write the final result
    if (tidx == 0) {
        v_r[tid / group_size] = shared[0];
    }
}

