//
//  curve_config.metal
//  Test
//
//  Created by Taras Shchybovyk on 4.11.2023.
//
#pragma once

#include "../../primitives/field.metal"
#include "../../primitives/projective.metal"

#include "params.metal"

namespace BN254 {
    typedef Field<PARAMS_BN254::fp_config> scalar_t;
    typedef Field<PARAMS_BN254::fp_config> point_field_t;
    static constant constexpr point_field_t gen_x = point_field_t{PARAMS_BN254::g1_gen_x};
    static constant constexpr point_field_t gen_y = point_field_t{PARAMS_BN254::g1_gen_y};
    static constant constexpr point_field_t b = point_field_t{PARAMS_BN254::weierstrass_b};
    typedef Projective<point_field_t, scalar_t, b, gen_x, gen_y> projective_t;
    typedef Affine<point_field_t> affine_t;
}
