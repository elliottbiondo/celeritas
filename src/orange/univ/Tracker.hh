//----------------------------------*-C++-*----------------------------------//
// Copyright 2021-2023 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file orange/univ/Tracker.hh
//---------------------------------------------------------------------------//
#pragma once

#include "corecel/Assert.hh"
#include "corecel/data/HyperslabIndexer.hh"
#include "corecel/grid/NonuniformGrid.hh"
#include "corecel/math/Algorithms.hh"
#include "orange/OrangeData.hh"
#include "orange/univ/detail/RaggedRightIndexer.hh"

#include "detail/Types.hh"
#include "detail/Utils.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Tracker base class
 */
class Tracker
{
  public:
    //!@{
    //! \name Type aliases
    using Initialization = detail::Initialization;
    using Intersection = detail::Intersection;
    using LocalState = detail::LocalState;

  public:
    CELER_FUNCTION Tracker(){};  // Default constructor

    //// ACCESSORS ////

    // DaughterId of universe embedded in a given volume
    virtual CELER_FUNCTION DaughterId daughter(LocalVolumeId vol) const = 0;

    //// OPERATIONS ////

    // Find the local volume from a position
    virtual CELER_FUNCTION Initialization initialize(LocalState const& state) const
        = 0;

    // Find the new volume by crossing a surface
    virtual CELER_FUNCTION Initialization
    cross_boundary(LocalState const& state) const
        = 0;

    // Calculate the distance to an exiting face for the current volume
    virtual CELER_FUNCTION Intersection intersect(LocalState const& state) const
        = 0;

    // Calculate nearby distance to an exiting face for the current volume
    virtual CELER_FUNCTION Intersection intersect(LocalState const& state,
                                                  real_type max_dist) const
        = 0;

    // Calculate closest distance to a surface in any direction
    virtual CELER_FUNCTION real_type safety(Real3 const& pos,
                                            LocalVolumeId vol) const
        = 0;

    // Calculate the local surface normal
    virtual CELER_FUNCTION Real3 normal(Real3 const& pos,
                                        LocalSurfaceId surf) const
        = 0;
};

//---------------------------------------------------------------------------//
}  // namespace celeritas
