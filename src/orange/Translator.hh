//----------------------------------*-C++-*----------------------------------//
// Copyright 2021-2022 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file orange/Translator.hh
//---------------------------------------------------------------------------//
#pragma once

#include "corecel/cont/Range.hh"
#include "orange/Types.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Class for performing spatial translations of 3D points
 */
class Translator
{
  public:
    // Constructor
    CELER_FUNCTION inline Translator(){};

    // Translate a single point
    CELER_FORCEINLINE_FUNCTION void
    operator()(Real3& point, const Real3& trans_vec) const;
};

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
/*!
 * Translate a single point
 */
CELER_FORCEINLINE_FUNCTION void
Translator::operator()(Real3& point, const Real3& trans_vec) const
{
    for (const auto& i : range(point.size()))
    {
        point[i] += trans_vec[i];
    }
}

//---------------------------------------------------------------------------//
} // namespace celeritas
