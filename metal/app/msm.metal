//
//  msm.metal
//  Test
//
//  Created by Taras Shchybovyk on 31.10.2023.
//
//
#include <metal_stdlib>
using namespace metal;


template <typename T>
kernel void ssm_kernel(device T* scalars [[buffer(0)]],
                       device T* points [[buffer(1)]],
                       device T* results [[buffer(2)]],
                       constant uint& N [[buffer(3)]],
                       uint tid [[thread_position_in_grid]]) {
    if (tid < N) {
        results[tid] = scalars[tid] * points[tid];
    }
}

template <typename P, typename A>
kernel void to_proj_kernel(device A* affine_points [[buffer(0)]],
                           device P* proj_points [[buffer(1)]],
                           constant uint& N [[buffer(3)]],
                           uint tid [[thread_position_in_grid]]) {
    if (tid < N) {
        proj_points[tid] = P::from_affine(affine_points[tid]);
    }
}

template <typename P>
kernel void sum_reduction_kernel(device P* v [[buffer(0)]],
                                 device P* v_r [[buffer(1)]],
                                 uint group_size [[threads_per_threadgroup]],
                                 uint tid [[thread_position_in_grid]],
                                 uint tidx [[thread_index_in_threadgroup]]) {
    // Perform first level of reduction:
    // Each thread loads one element from global to shared mem
    threadgroup P shared[1024];
    
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
