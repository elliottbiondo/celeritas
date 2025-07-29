//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file orange/orangeinp/detail/SpecialTrapezoidDecomposer.hh
//---------------------------------------------------------------------------//
#pragma once

#include <algorithm>

#include "corecel/Assert.hh"
#include "corecel/Types.hh"
#include "corecel/cont/Array.hh"
#include "corecel/cont/Range.hh"
#include "corecel/math/SoftEqual.hh"
#include "orange/OrangeTypes.hh"
#include "orange/orangeinp/CsgTypes.hh"
#include "orange/univ/detail/Utils.hh"

namespace
{
//! Convenience enumeration for implementations in this file
enum
{
    R = 0,
    Z = 1
};

constexpr auto left = Bound::lo;
constexpr auto right = Bound::hi;

}  // namespace

namespace celeritas
{
namespace orangeinp
{
namespace detail
{
//---------------------------------------------------------------------------//
/*!
 * Decompose a polygon into a
 *
 */
class SpecialTrapezoidDecomposer
{
  public:
    //!@{
    //! \name Type aliases
    using Real2 = celeritas::Array<real_type, 2>;
    using VecReal2 = std::vector<Real2>;
    using VecSpecTrap = std::vector<SpecialTrapezoid>;
    //!@}

  public:
    // Construct with vector of ordered points and a tolerance
    explicit SpecialTrapezoidDecomposer(VecReal2 points,
                                        Tolerance<> const& tol);

    // Decompose polygon into special trapezoids
    VecSpecTrap operator()();

  private:
    /// TYPES ///
    using InclusiveRange = EnumArray<Bound, size_type>;

    celeritas::Array<size_type, 2>;

    /// DATA ///
    VecReal2 points_;
    SoftEqual<real_type> soft_equal_;

    /// HELPER FUNCTIONS ///

    // Clip off the special trapezoid created
    SpecialTrapezoid clip_lowest_trap();

    // Delete all points in the range [left, right].
    void delete_range(InclusiveRange inclusive_range);

    // Create a SoftEqual object for internal use
    SoftEqual<real_type>
    make_soft_equal(VecReal2 polygon, Tolerance<> const& tol) const;

    // Determine the next index, with modular indexing
    size_type calc_next(size_type i) const;

    // Determine the previous index, with modular indexing
    size_type calc_previous(size_type i) const;
};

//---------------------------------------------------------------------------//
/*!
 * Construct with vector of ordered points.
 */
SpecialTrapezoidDecomposer::SpecialTrapezoidDecomposer(VecReal2 points,
                                                       Tolerance<> const& tol)
    : points_{points}, soft_equal_(this->make_soft_equal(points, tol))
{
    CELER_EXPECT(points_.size() > 2);
}

//---------------------------------------------------------------------------//
/*!
 * Decompose polygon into special trapezoids
 */
auto SpecialTrapezoidDecomposer::operator()() -> VecSpecTrap
{
    VecSpecTrap result;

    while (points_.size() > 2)
    {
        result.push_back(this->clip_lowest_trap());
    }

    return result;
}

//---------------------------------------------------------------------------//
/*!
 * Clip off the special trapezoid created.
 *
 */
SpecialTrapezoid SpecialTrapezoidDecomposer::clip_lowest_trap()
{
    auto const& p = points_;

    auto br = this->calc_bottom_range();

    SpecialTrapezoid::ZSegment zseg_bot;
    zseg_bot.r = {p[br[left]][R], p[br[right]][R]};
    zeg_bot.z = p[br[left]][Z];

    // Upper left and upper right candidates for the minimum z
    ul = this->calc_prev(br[left]);
    ur = this->calc_next(br[right]);

    SpecialTrapezoid::ZSegment zseg_top;
    if (soft_equal_(p[ul][Z], p[ur][Z]))
    {
        zeg_top.r = {p[ul][R], p[ur][R]};
        zeg_top.z = p[ul][Z];
    }
    else if (p[ul][Z] < p[ur][Z])
    {
        auto new_r = this->interpolate();
        zeg_top.r = {p[ul][R], p[ur][R]};
        zeg_top.z = p[ul][Z];
    }

    this->erase_range();
    return {zseg_bot, zseg_top};
}

}
//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * Find the r value on the line between p0 and p1 with a supplied z value.
 */
real_type SpecialTrapezoidDecomposer::interpolate(Real2 const& p0,
                                                  Real2 const& p1,
                                                  real_type z) const
{
    return p0[R] + (z - p0[Z]) / (p1[Z] - p0[Z]) * (p1[R] - p0[R]);
}

//---------------------------------------------------------------------------//
/*!
 * Find a contiguous range of values that all have the minimum z value.
 */
InclusiveRange SpecialTrapezoidDecomposer::calc_bottom_range() const
{
    // Find the index of the point with the lowest z
    auto starting_it = std::min_element(
        points_.begin(), points_.end(), [](Real2 const& a, Real2 const& b) {
            return a[1] < b[1];
        });
    auto min_idx = std::distance(points_.begin(), starting_it);

    // Traverse the points by applying the step function until a point is
    // encountered that does not have the minimum z
    auto traverse_min_z = [&](auto step) {
        size_type i = min_idx;
        while (soft_equal_(points_[min_idx][Z], points_[step(i)][Z]))
        {
            i = step(i);
        }
        return i;
    };

    // Find the leftmost and rightmost points that have the minimum z
    size_type bot_left = traverse_min_z(this->calc_previous);
    size_type bot_right = traverse_min_z(this->calc_next);

    return {bot_left, bot_right};
}

//---------------------------------------------------------------------------//
/*!
 * Delete all points in the inclusive range of indices.
 */
void SpecialTrapezoidDecomposer::delete_range(InclusiveRange ir)
{
    points_.erase(points_.begin() + ir[left], points_.begin() + ir[right] + 1);
}

//---------------------------------------------------------------------------//
/*!
 * Make a SoftEqual object based on polygon extents and a tolerance
 */
SoftEqual<real_type>
ConvexHullFinder<T>::make_soft_equal(VecReal2 const& polygon,
                                     Tolerance<> const& t) const
{
    auto const [r_min, r_max]
        = find_extrema(polygon, static_cast<size_type>(R));
    auto const [z_min, z_max]
        = find_extrema(polygon, static_cast<size_type>(Z));

    // Convert min/max x and y values to extents
    Real3 const extents{r_max - r_min, z_max - z_min, 0};
    return SoftOrientation<real_type>(
        ::celeritas::detail::BumpCalculator(t)(extents));
}

//---------------------------------------------------------------------------//
/*!
 * Determine the next index using modular indexing.
 */
size_type SpecialTrapezoidDecomposer::calc_next(size_type i) const
{
    return (i + 1) % points_.size();
}

//---------------------------------------------------------------------------//
/*!
 * Determine the previous index using modular indexing.
 */
size_type SpecialTrapezoidDecomposer::calc_previous(size_type i) const
{
    return (i - 1) % points_.size();
}

//---------------------------------------------------------------------------//
}  // namespace detail
}  // namespace orangeinp
}  // namespace celeritas
