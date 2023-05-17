//----------------------------------*-C++-*----------------------------------//
// Copyright 2021-2023 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file orange/univ/RectArrayTracker.hh
//---------------------------------------------------------------------------//
#pragma once

#include "corecel/Assert.hh"
#include "corecel/data/HyperslabIndexer.hh"
#include "corecel/data/RaggedRightIndexer.hh"
#include "corecel/grid/NonuniformGrid.hh"
#include "corecel/math/Algorithms.hh"
#include "orange/OrangeData.hh"

#include "detail/Types.hh"
#include "detail/Utils.hh"

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Track a particle within an axes-aligned rectilinear grid.
 */
class RectArrayTracker
{
  public:
    //!@{
    //! \name Type aliases
    using ParamsRef = NativeCRef<OrangeParamsData>;
    using Initialization = detail::Initialization;
    using Intersection = detail::Intersection;
    using LocalState = detail::LocalState;
    using Grid = NonuniformGrid<real_type>;
    using VolumeIndexer = HyperslabIndexer<3>;
    using VolumeInverseIndexer = HyperslabInverseIndexer<3>;
    using SurfaceIndexerData = RaggedRightIndexerData<3>;
    using SurfaceIndexer = RaggedRightIndexer<3>;
    using SurfaceInverseIndexer = RaggedRightInverseIndexer<3>;
    using Coords = Array<size_type, 3>;
    //!@}

  public:
    // Construct with parameters (unit definitions and this one's ID)
    inline CELER_FUNCTION
    RectArrayTracker(ParamsRef const& params, RectArrayId rid);

    //// ACCESSORS ////

    //! Number of local volumes
    CELER_FUNCTION LocalVolumeId::size_type num_volumes() const
    {
        return record_.daughters.size();
    }

    //! Number of local surfaces
    CELER_FUNCTION LocalSurfaceId::size_type num_surfaces() const
    {
        size_type num_surfs = 0;
        for (auto ax : range(Axis::size_))
        {
            num_surfs += dims_[to_int(ax)] + 1;
        }
        return num_surfs;
    }

    //! RectArrayRecord for this tracker
    CELER_FUNCTION RectArrayRecord const& record() const { return record_; }

    // DaughterId of universe embedded in a given volume
    inline CELER_FUNCTION DaughterId daughter(LocalVolumeId vol) const;

    ////// OPERATIONS ////

    // Find the local volume from a position
    inline CELER_FUNCTION Initialization
    initialize(LocalState const& state) const;

    // Calculate distance-to-intercept for the next surface
    inline CELER_FUNCTION Intersection intersect(LocalState const& state) const;

    // Calculate distance-to-intercept for the next surface, with max distance
    inline CELER_FUNCTION Intersection intersect(LocalState const& state,
                                                 real_type max_dist) const;

    // Find the new volume by crossing a surface
    inline CELER_FUNCTION Initialization
    cross_boundary(LocalState const& state) const;

    // Calculate closest distance to a surface in any direction
    inline CELER_FUNCTION real_type safety(Real3 const& pos,
                                           LocalVolumeId vol) const;

    // Calculate the local surface normal
    inline CELER_FUNCTION Real3 normal(Real3 const& pos,
                                       LocalSurfaceId surf) const;

  private:
    //// DATA ////
    ParamsRef const& params_;
    RectArrayRecord const& record_;
    Array<size_type, 3> dims_;
    SurfaceIndexerData surface_indexer_data_;

    //// METHODS ////

    template<class F>
    inline CELER_FUNCTION Intersection intersect_impl(LocalState const&,
                                                      F) const;
};

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
/*!
 * Construct with reference to persistent parameter data.
 */
CELER_FUNCTION
RectArrayTracker::RectArrayTracker(ParamsRef const& params, RectArrayId rid)
    : params_(params), record_(params.rect_arrays[rid])
{
    CELER_EXPECT(params_);

    for (auto ax : range(Axis::size_))
    {
        dims_[to_int(ax)] = record_.grid[to_int(ax)].size() - 1;
    }

    SurfaceIndexerData::Sizes sizes{dims_[to_int(Axis::x)] + 1,
                                    dims_[to_int(Axis::y)] + 1,
                                    dims_[to_int(Axis::z)] + 1};

    surface_indexer_data_ = SurfaceIndexerData(sizes);
}

//---------------------------------------------------------------------------//
/*!
 * Find the local volume from a position.
 *
 * To avoid edge cases and inconsistent logical/physical states, it is
 * prohibited to initialize from an arbitrary point directly onto a surface.
 */
CELER_FUNCTION auto RectArrayTracker::initialize(LocalState const& state) const
    -> Initialization
{
    CELER_EXPECT(params_);
    CELER_EXPECT(!state.surface && !state.volume);

    Coords coords;

    for (auto ax : range(Axis::size_))
    {
        auto const& pos = state.pos[to_int(ax)];
        Grid grid(record_.grid[to_int(ax)], params_.reals);

        if (pos < grid.front() || pos > grid.back())
        {
            // Outside the rect array
            return {record_.background, {}};
        }
        else
        {
            size_type index = grid.find(pos);
            bool edge = grid[index] == pos;

            if (!edge)
            {
                coords[to_int(ax)] = index;
            }
            else
            {
                // On boundary
                return {record_.background, {}};
            }
        }
    }

    VolumeIndexer vi(dims_);
    return {LocalVolumeId{vi(coords)}, {}};
}

//---------------------------------------------------------------------------//
/*!
 * Find the local volume on the opposite side of a surface.
 */
CELER_FUNCTION auto
RectArrayTracker::cross_boundary(LocalState const& state) const
    -> Initialization
{
    CELER_EXPECT(state.surface && state.volume);

    // Find the coords of the current volume
    VolumeIndexer vi(dims_);
    VolumeInverseIndexer vii(dims_);
    auto coords = vii(state.volume.unchecked_get());

    // Find the index of axis (x/y/z) we are about to cross:
    SurfaceInverseIndexer sii(surface_indexer_data_);
    auto ax_idx = sii(state.surface.id().unchecked_get())[0];

    // Value for incrementing the axial coordinate upon crossing
    int inc = (state.surface.sense() == Sense::outside) ? -1 : 1;

    detail::OnLocalSurface new_surface(state.surface.id(),
                                       flip_sense(state.surface.sense()));

    if ((coords[ax_idx] == 0 && inc == -1)
        || (coords[ax_idx] == dims_[ax_idx] - 1 && inc == 1))
    {
        // crossimg out
        return {record_.background, new_surface};
    }
    else
    {
        coords[ax_idx] += inc;
        return {LocalVolumeId(vi(coords)), new_surface};
    }
}

//---------------------------------------------------------------------------//
/*!
 * Calculate distance-to-intercept for the next surface.
 */
CELER_FUNCTION auto RectArrayTracker::intersect(LocalState const& state) const
    -> Intersection
{
    Intersection result = this->intersect_impl(state, detail::IsFinite{});
    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Calculate distance-to-intercept for the next surface, with max distance.
 */
CELER_FUNCTION auto
RectArrayTracker::intersect(LocalState const& state, real_type max_dist) const
    -> Intersection
{
    CELER_EXPECT(max_dist > 0);
    Intersection result
        = this->intersect_impl(state, detail::IsNotFurtherThan{max_dist});
    if (!result)
    {
        result.distance = max_dist;
    }
    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Calculate nearest distance to a surface in any direction.
 *
 * On an axis-aligned rectlinear grid the minimum distance to any surface is
 * always always occurs along a line parallel to an axis.
 */
CELER_FUNCTION real_type RectArrayTracker::safety(Real3 const& pos,
                                                  LocalVolumeId volid) const
{
    CELER_EXPECT(volid && volid.get() < this->num_volumes());

    VolumeInverseIndexer vii(dims_);
    auto coords = vii(volid.unchecked_get());

    real_type min_dist = numeric_limits<real_type>::infinity();

    for (auto ax : range(Axis::size_))
    {
        for (auto i : range(2))
        {
            auto target_coord = coords[to_int(ax)] + i;
            real_type target
                = params_.reals[record_.grid[to_int(ax)]][target_coord];
            min_dist = min(min_dist, fabs(pos[to_int(ax)] - target));
        }
    }

    CELER_ENSURE(min_dist >= 0
                 && min_dist < numeric_limits<real_type>::infinity());
    return min_dist;
}

//---------------------------------------------------------------------------//
/*!
 * Calculate the local surface normal.
 */
CELER_FUNCTION auto RectArrayTracker::normal([[maybe_unused]] Real3 const& pos,
                                             LocalSurfaceId surf) const -> Real3
{
    CELER_EXPECT(surf && surf.get() < this->num_surfaces());
    SurfaceInverseIndexer sii(surface_indexer_data_);
    size_type ax = sii(surf.get())[0];

    Real3 normal{0., 0., 0.};
    normal[ax] = 1.0;

    return normal;
}

//---------------------------------------------------------------------------//
/*!
 * DaughterId of universe embedded in a given volume.
 */
CELER_FORCEINLINE_FUNCTION DaughterId
RectArrayTracker::daughter(LocalVolumeId vol) const
{
    CELER_EXPECT(vol && vol.get() < this->num_volumes());
    return record_.daughters[vol];
}

//---------------------------------------------------------------------------//
// PRIVATE INLINE DEFINITIONS
//---------------------------------------------------------------------------//
/*!
 * Calculate distance-to-intercept for the next surface.
 */
template<class F>
CELER_FUNCTION auto
RectArrayTracker::intersect_impl(LocalState const& state, F is_valid) const
    -> Intersection
{
    CELER_EXPECT(state.volume && !state.temp_sense.empty());

    VolumeInverseIndexer vii(dims_);
    auto coords = vii(state.volume.unchecked_get());

    Intersection result;
    Sense sense;
    SurfaceIndexer si(surface_indexer_data_);

    for (auto ax : range(Axis::size_))
    {
        auto dir = state.dir[to_int(ax)];

        // Ignore any stationary axis
        if (dir == 0.0)
        {
            continue;
        }

        auto target_coord = coords[to_int(ax)] + static_cast<int>(dir > 0.);

        auto target_value
            = params_.reals[record_.grid[to_int(ax)]][target_coord];

        double dist
            = (target_value - static_cast<double>(state.pos[to_int(ax)]))
              / state.dir[to_int(ax)];

        if (dist > 0 && dist < result.distance)
        {
            result.distance = dist;

            sense = dir > 0 ? Sense::inside : Sense::outside;
            auto local_surface = LocalSurfaceId(
                si({static_cast<size_type>(to_int(ax)), target_coord}));
            result.surface = detail::OnLocalSurface(local_surface, sense);
        }
    }

    return (is_valid(result.distance)) ? result : Intersection{};
}

//---------------------------------------------------------------------------//
}  // namespace celeritas
