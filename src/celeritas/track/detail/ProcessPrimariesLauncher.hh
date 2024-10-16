//----------------------------------*-C++-*----------------------------------//
// Copyright 2020-2022 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celeritas/track/detail/ProcessPrimariesLauncher.hh
//---------------------------------------------------------------------------//
#pragma once

#include "corecel/cont/Span.hh"

#include "../TrackInitData.hh"

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
/*!
 * Create track initializers from primary particles.
 */
template<MemSpace M>
class ProcessPrimariesLauncher
{
  public:
    //!@{
    //! Type aliases
    using TrackInitStateRef = TrackInitStateData<Ownership::reference, M>;
    //!@}

  public:
    // Construct with shared and state data
    CELER_FUNCTION ProcessPrimariesLauncher(Span<const Primary>      primaries,
                                            const TrackInitStateRef& data)
        : primaries_(primaries), data_(data)
    {
        CELER_EXPECT(data_);
    }

    // Create track initializers from primaries
    inline CELER_FUNCTION void operator()(ThreadId tid) const;

  private:
    Span<const Primary>      primaries_;
    const TrackInitStateRef& data_;
};

//---------------------------------------------------------------------------//
/*!
 * Create track initializers from primaries.
 */
template<MemSpace M>
CELER_FUNCTION void ProcessPrimariesLauncher<M>::operator()(ThreadId tid) const
{
    TrackInitializer& init    = data_.initializers[ThreadId(
        data_.initializers.size() - primaries_.size() + tid.get())];
    const Primary&    primary = primaries_[tid.get()];

    // Construct a track initializer from a primary particle
    init.sim.track_id         = primary.track_id;
    init.sim.parent_id        = TrackId{};
    init.sim.event_id         = primary.event_id;
    init.sim.num_steps        = 0;
    init.sim.time             = primary.time;
    init.sim.status           = TrackStatus::alive;
    init.geo.pos              = primary.position;
    init.geo.dir              = primary.direction;
    init.particle.particle_id = primary.particle_id;
    init.particle.energy      = primary.energy;
}

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas
