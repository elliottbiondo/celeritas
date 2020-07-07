//---------------------------------*-CUDA-*----------------------------------//
// Copyright 2020 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file RngEngine.i.cuh
//---------------------------------------------------------------------------//

#include <cstddef>

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Construct from state
 */
__device__ RngEngine::RngEngine(RngState* state) : state_(state) {}

//---------------------------------------------------------------------------//
/*!
 * Sample a random number
 */
__device__ auto RngEngine::operator()() -> result_type
{
    return curand(state_);
}

//---------------------------------------------------------------------------//
} // namespace celeritas
