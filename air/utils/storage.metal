//
//  storage.metal
//  Test
//
//  Created by Taras Shchybovyk on 31.10.2023.
//

#include <metal_stdlib>
using namespace metal;


#define LIMBS_ALIGNMENT(x) ((x) % 4 == 0 ? 16 : ((x) % 2 == 0 ? 8 : 4))

template <unsigned LIMBS_COUNT>
struct storage
{
  static constant constexpr unsigned LC = LIMBS_COUNT;
  alignas(LIMBS_ALIGNMENT(LIMBS_COUNT)) uint32_t limbs[LIMBS_COUNT];
};

template <unsigned OMEGAS_COUNT, unsigned LIMBS_COUNT>
struct storage_array
{
  alignas(LIMBS_ALIGNMENT(LIMBS_COUNT)) storage<LIMBS_COUNT> storages[OMEGAS_COUNT];
};