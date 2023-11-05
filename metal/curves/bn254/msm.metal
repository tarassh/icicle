//
//  msm.metal
//  Test
//
//  Created by Taras Shchybovyk on 5.11.2023.
//

#pragma once

#include "curve_config.metal"
#include "../../app/msm.metal"
#include <metal_stdlib>
using namespace metal;

kernel
void to_proj_kernel_bn254(device BN254::projective_t *proj_points [[buffer(0)]],
                          device BN254::affine_t *affine_points [[buffer(1)]],
                          device const uint &N [[buffer(2)]],
                          uint tid [[thread_position_in_grid]])
{
    
    to_proj_kernel(affine_points, proj_points, N, tid);
}

kernel
void ssm_kernel_bn254(device BN254::scalar_t *scalars [[buffer(0)]],
                      device BN254::projective_t *points [[buffer(1)]],
                      device BN254::projective_t *results [[buffer(2)]],
                      device const uint &N [[buffer(3)]],
                      uint tid [[thread_position_in_grid]])
{
    ssm_kernel(scalars, points, results, N, tid);
}

kernel
void sum_reduction_kernel_bn254(device BN254::projective_t* v [[buffer(0)]],
                                device BN254::projective_t* v_r [[buffer(1)]],
                                uint group_size [[threads_per_threadgroup]],
                                uint tid [[thread_position_in_grid]],
                                uint tidx [[thread_index_in_threadgroup]]) {
    // TODO: proper size for shared
    threadgroup BN254::projective_t shared[1024];
    
    sum_reduction_kernel(v, v_r, shared, group_size, tid, tidx);
}
