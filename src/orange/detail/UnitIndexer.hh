//----------------------------------*-C++-*----------------------------------//
// Copyright 2022 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file orange/detail/UnitIndexer.hh
//---------------------------------------------------------------------------//
#pragma once

#include <vector>

#include "orange/Types.hh"

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
/*!
 * Convert a unit input to params data.
 *
 * Linearize the data in a UnitInput and add it to the host.
 */
class UnitIndexer
{
  public:
    //!@{
    //! \name Type aliases
    using Local_Surface = std::tuple<UniverseId, SurfaceId>;
    using Local_Volume  = std::tuple<UniverseId, VolumeId>;
    using Vec_Size      = std::vector<size_type>;
    //!@}

  public:
    // Construct from sizes
    UnitIndexer(Vec_Size num_surfaces, Vec_Size num_volumes);

    // Local-to-global
    inline SurfaceId global_surface(UniverseId uni, SurfaceId surface) const;
    inline VolumeId  global_volume(UniverseId uni, VolumeId volume) const;

    // Global-to-local
    inline Local_Surface local_surface(SurfaceId id) const;
    inline Local_Volume  local_volume(VolumeId id) const;

    //! Total number of universes
    size_type num_universes() const { return d_surfaces.size() - 1; }

    //! Total number of surfaces
    size_type num_surfaces() const { return d_surfaces.back(); }

    //! Total number of cells
    size_type num_volumes() const { return d_volumes.back(); }

  private:
    // >>> DATA
    Vec_Size d_surfaces;
    Vec_Size d_volumes;

    // >>> IMPLEMENTATION METHODS
    static inline Vec_Size::const_iterator
    find_local(const Vec_Size& offsets, size_type id);

    static inline size_type local_size(const Vec_Size& offsets, UniverseId uni);
};

//---------------------------------------------------------------------------//
} // namespace detail
} // namespace celeritas
