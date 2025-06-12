//------------------------------- -*- C++ -*- -------------------------------//
// Copyright Celeritas contributors: see top-level COPYRIGHT file for details
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file orange/orangeinp/detail/PolylineUtils.hh
//! \brief Utility functions for polylines in 2D or 3D space.
//---------------------------------------------------------------------------//
#pragma once

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "corecel/cont/Range.hh"

namespace celeritas
{
namespace orangeinp
{
namespace detail
{

//---------------------------------------------------------------------------//
/*
 * Check if polyline is monotonically increasing in the given dimension.
 *
 * Here, "strict" means that equal values are considered non-monotonoic
 */
template<class T, std::size_t N>
bool is_monotonic_increasing(std::vector<T> const& points,
                             std::size_t dim,
                             bool strict = true)
{
    CELER_EXPECT(points.size() > 1);
    CELER_EXPECT(dim < N);
    for (auto i : range<std::size_t>(points.size() - 1))
    {
        T curr = points[i][dim];
        T next = points[i + 1][dim];
        if ((strict && next <= curr) || (!strict && next < curr))
        {
            return false;
        }
    }
    return true;
}

//---------------------------------------------------------------------------//
/*
 * Decompose a polyline into its monotonic and non-monotonic subsequences.
 *
 * In the example below, the monotonic sequence is A-B-C-E-F-H and the
 * non-monotonic sequences are C-D-E and F-G-H.
 * \verbatim
                   . H
                 .  |
               .    |
      E  ____ F     |
        | .    \    |
        | .     \   |
        |  . C    \ |
        |  /\      \| G
        |/    \ B
        D      |
               |
               |
               A
   \endverbatim
 */
template<class T, size N>
class MonotonicDecomposer
{
  public:
    /// TYPES
    using VecT = std::vector<T>;

    struct ResultType
    {
        VecT monotonic_pl;
        std::vector<VecT> non_monotonic_pls;
    }

    //! Constructor
    ResultType
    MonotonicDecomposer(VecT const& points, size_type dim)
        : points_(points), dim(d)
    {
        CELER_EXPECT(dim < N);
        CELER_EXPECT(points.size() > 1);
    }

    //! Identify constituant monotonic points and non-monotonic polylines
    ResultType operator()() const
    {
        Result_type results;
        results.montonic_pls.append(points[0]);

        size_type i = 1 while (i < points_.size())
        {
            if (points[i][dim] < points[i - 1][dim])
            {
                results.monotonic_pls.append(points[i]);
            }
            else
            {
                results.monotonic_pls.append(points[i - 1]);
                do
                {
                    results.monotonic_pls.append(points[i]);
                }
            }
        }

        return results;
    }

  private:
    /// DATA ///
    VecT const& points;
    size_type dim;
};

//---------------------------------------------------------------------------//
}  // namespace detail
}  // namespace orangeinp
}  // namespace celeritas
