//----------------------------------*-C++-*----------------------------------//
// Copyright 2022-2024 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file orange/detail/BIHTraverser.hh
//---------------------------------------------------------------------------//
#pragma once

#include "corecel/math/Algorithms.hh"

#include "../BoundingBoxUtils.hh"
#include "../OrangeData.hh"

namespace celeritas
{
namespace detail
{
//---------------------------------------------------------------------------//
/*!
 * Traverse BIH tree using a depth-first search.
 *
 * \todo move to top-level orange directory out of detail namespace
 */
class BIHTraverser
{
  public:
    //!@{
    //! \name Type aliases
    using Storage = NativeCRef<BIHTreeData>;
    //!@}

    // Construct from vector of bounding boxes and storage for LocalVolumeIds
    inline CELER_FUNCTION
    BIHTraverser(BIHTree const& tree, Storage const& storage);

    // Traverse the tree
    template<class F>
    inline CELER_FUNCTION LocalVolumeId
    find_enclosing_volume(Real3 const& point, F&& visit_vol) const;

  private:
    //// TYPES ////

    //// DATA ////
    BIHTree const& tree_;
    Storage const& storage_;
    size_type leaf_offset_;

    //// HELPER FUNCTIONS ////

    // Traverse the tree
    template<class F, class G, class H>
    inline CELER_FUNCTION LocalVolumeId traverse(Real3 const& point,
                                                 F&& visit_edge,
                                                 G&& visit_bbox,
                                                 H&& visit_vol) const;

    // Get the ID of the next node in the traversal sequence
    template<class F>
    inline CELER_FUNCTION BIHNodeId next_node(BIHNodeId const& current_id,
                                              BIHNodeId const& previous_id,
                                              Real3 const& point,
                                              F&& visit_edge) const;

    // Determine if a node is inner, i.e., not a leaf
    inline CELER_FUNCTION bool is_inner(BIHNodeId id) const;

    // Get an inner node for a given BIHNodeId
    inline CELER_FUNCTION BIHInnerNode const&
    get_inner_node(BIHNodeId id) const;

    // Get a leaf node for a given BIHNodeId
    inline CELER_FUNCTION BIHLeafNode const& get_leaf_node(BIHNodeId id) const;

    // Determine if any leaf node volumes contain the point
    template<class F, class G>
    inline CELER_FUNCTION LocalVolumeId visit_leaf(BIHLeafNode const& leaf_node,
                                                   Real3 const& point,
                                                   F&& visit_bbox,
                                                   G&& visit_vol) const;

    // Determine if any inf_vols contain the point
    template<class F>
    inline CELER_FUNCTION LocalVolumeId visit_inf_vols(F&& visit_vol) const;

    // Determine if a single bbox contains the point
    inline CELER_FUNCTION bool
    visit_bbox(LocalVolumeId const& id, Real3 const& point) const;
};

//---------------------------------------------------------------------------//
// INLINE DEFINITIONS
//---------------------------------------------------------------------------//
/*!
 * Construct from vector of bounding boxes and storage.
 */
CELER_FUNCTION
BIHTraverser::BIHTraverser(BIHTree const& tree,
                           BIHTraverser::Storage const& storage)
    : tree_(tree), storage_(storage), leaf_offset_(tree.inner_nodes.size())
{
    CELER_EXPECT(tree);
}

//---------------------------------------------------------------------------//
/*!
 * Point-in-volume operation.
 */
template<class F>
CELER_FUNCTION LocalVolumeId
BIHTraverser::find_enclosing_volume(Real3 const& point, F&& visit_vol) const
{
    auto visit_edge = [](BIHInnerNode const& node,
                         BIHInnerNode::Side side,
                         Real3 const& point) {
        CELER_EXPECT(side < BIHInnerNode::Side::size_);

        auto pos = node.edges[side].bounding_plane_pos;
        auto point_pos = point[to_int(node.axis)];

        return (side == BIHInnerNode::Side::left) ? (point_pos < pos)
                                                  : (pos < point_pos);
    };

    auto visit_bbox = [this](LocalVolumeId const& id, Real3 const& point) {
        return is_inside(storage_.bboxes[tree_.bboxes[id]], point);
    };

    return this->traverse(point, visit_edge, visit_bbox, visit_vol);
}

//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * Point-in-volume operation.
 */
template<class F, class G, class H>
CELER_FUNCTION LocalVolumeId BIHTraverser::traverse(Real3 const& point,
                                                    F&& visit_edge,
                                                    G&& visit_bbox,
                                                    H&& visit_vol) const
{
    BIHNodeId previous_node;
    BIHNodeId current_node{0};
    LocalVolumeId id;

    do
    {
        if (!this->is_inner(current_node))
        {
            id = this->visit_leaf(
                this->get_leaf_node(current_node), point, visit_bbox, visit_vol);

            if (id)
            {
                return id;
            }
        }

        previous_node = exchange(
            current_node,
            this->next_node(current_node, previous_node, point, visit_edge));

    } while (current_node);

    if (!id)
    {
        id = this->visit_inf_vols(visit_vol);
    }

    return id;
}

//---------------------------------------------------------------------------//
/*!
 *  Get the ID of the next node in the traversal sequence.
 */
template<class F>
CELER_FUNCTION BIHNodeId BIHTraverser::next_node(BIHNodeId const& current_id,
                                                 BIHNodeId const& previous_id,
                                                 Real3 const& point,
                                                 F&& visit_edge) const
{
    using Side = BIHInnerNode::Side;

    BIHNodeId next_id;

    if (this->is_inner(current_id))
    {
        auto const& current_node = this->get_inner_node(current_id);
        if (previous_id == current_node.parent)
        {
            // Visiting this inner node for the first time; go down either left
            // or right edge
            if (visit_edge(current_node, Side::left, point))
            {
                next_id = current_node.edges[Side::left].child;
            }
            else
            {
                next_id = current_node.edges[Side::right].child;
            }
        }
        else if (previous_id == current_node.edges[Side::left].child)
        {
            // Visiting this inner node for the second time; go down right edge
            // or return to parent
            if (visit_edge(current_node, Side::right, point))
            {
                next_id = current_node.edges[Side::right].child;
            }
            else
            {
                next_id = current_node.parent;
            }
        }
        else
        {
            // Visiting this inner node for the third time; return to parent
            CELER_EXPECT(previous_id == current_node.edges[Side::right].child);
            next_id = current_node.parent;
        }
    }
    else
    {
        // Leaf node; return to parent
        CELER_EXPECT(previous_id == this->get_leaf_node(current_id).parent);
        next_id = previous_id;
    }

    return next_id;
}

//---------------------------------------------------------------------------//
/*!
 *  Determine if a node is inner, i.e., not a leaf.
 */
CELER_FUNCTION
bool BIHTraverser::is_inner(BIHNodeId id) const
{
    return id.unchecked_get() < leaf_offset_;
}

//---------------------------------------------------------------------------//
/*!
 *  Get an inner node for a given BIHNodeId.
 */
CELER_FUNCTION
BIHInnerNode const& BIHTraverser::get_inner_node(BIHNodeId id) const
{
    CELER_EXPECT(this->is_inner(id));
    return storage_.inner_nodes[tree_.inner_nodes[id.unchecked_get()]];
}

//---------------------------------------------------------------------------//
/*!
 *  Get a leaf node for a given BIHNodeId.
 */
CELER_FUNCTION
BIHLeafNode const& BIHTraverser::get_leaf_node(BIHNodeId id) const
{
    CELER_EXPECT(!this->is_inner(id));
    return storage_
        .leaf_nodes[tree_.leaf_nodes[id.unchecked_get() - leaf_offset_]];
}

//---------------------------------------------------------------------------//
/*!
 * Determine if any leaf node volumes contain the point.
 */
template<class F, class G>
CELER_FUNCTION LocalVolumeId
BIHTraverser::visit_leaf(BIHLeafNode const& leaf_node,
                         Real3 const& point,
                         F&& visit_bbox,
                         G&& visit_vol) const
{
    for (auto i : range(leaf_node.vol_ids.size()))
    {
        auto id = storage_.local_volume_ids[leaf_node.vol_ids[i]];
        if (visit_bbox(id, point) && visit_vol(id))
        {
            return id;
        }
    }
    return LocalVolumeId{};
}

//---------------------------------------------------------------------------//
/*!
 * Determine if any volumes in inf_vols contain the point.
 */
template<class F>
CELER_FUNCTION LocalVolumeId BIHTraverser::visit_inf_vols(F&& visit_vol) const
{
    for (auto i : range(tree_.inf_volids.size()))
    {
        auto id = storage_.local_volume_ids[tree_.inf_volids[i]];
        if (visit_vol(id))
        {
            return id;
        }
    }
    return LocalVolumeId{};
}

//---------------------------------------------------------------------------//
}  // namespace detail
}  // namespace celeritas
