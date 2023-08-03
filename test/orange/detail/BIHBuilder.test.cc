//----------------------------------*-C++-*----------------------------------//
// Copyright 2021-2023 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file orange/detail/BIHBuilder.test.cc
//---------------------------------------------------------------------------//
#include "orange/detail/BIHBuilder.hh"

#include "corecel/data/CollectionBuilder.hh"
#include "corecel/data/CollectionMirror.hh"
#include "orange/detail/BIHData.hh"
#include "celeritas/Types.hh"

#include "celeritas_test.hh"

using BIHBuilder = celeritas::detail::BIHBuilder;
using BIHInnerNode = celeritas::detail::BIHInnerNode;
using BIHLeafNode = celeritas::detail::BIHLeafNode;

namespace celeritas
{
namespace test
{
class BIHBuilderTest : public Test
{
  public:
    void SetUp()
    {
        auto inf = std::numeric_limits<fast_real_type>::infinity();
        bboxes_.push_back({{-inf, -inf, -inf}, {inf, inf, inf}});
    }

  protected:
    std::vector<FastBBox> bboxes_;

    BIHTreeData<Ownership::value, MemSpace::host> storage_;
};

//---------------------------------------------------------------------------//
/* Simple test with partial and fully overlapping bounding boxes.
 *
 *         0    V1    1.6
 *         |--------------|
 *
 *                    1.2   V2    2.8
 *                    |---------------|
 *    y=1 ____________________________________________________
 *        |           |   |           |                      |
 *        |           |   |           |         V3           |
 *    y=0 |___________|___|___________|______________________|
 *        |                                                  |
 *        |             V4, V5 (total overlap)               |
 *   y=-1 |__________________________________________________|
 *
 *        x=0                                                x=5
 *
 * Resultant tree structure in terms of BIHNodeIds (N) and volumes (V):
 *
 *                      ___ N0 ___
 *                    /            \
 *                  N1              N2
 *                 /  \           /    \
 *                N3   N4        N5     N6
 *                V1   V2       V4,V5   V3
 *
 * In terms of BIHInnerNodeIds (I) and BIHLeafNodeIds (L):
 *
 *                      ___ I0 ___
 *                    /            \
 *                  I1              I2
 *                 /  \           /    \
 *                L1   L2        L3     L4
 *                V1   V2       V4,V5   V3
 */
TEST_F(BIHBuilderTest, basic)
{
    using Edge = BIHInnerNode::Edge;

    bboxes_.push_back({{0, 0, 0}, {1.6, 1, 100}});
    bboxes_.push_back({{1.2, 0, 0}, {2.8, 1, 100}});
    bboxes_.push_back({{2.8, 0, 0}, {5, 1, 100}});
    bboxes_.push_back({{0, -1, 0}, {5, 0, 100}});
    bboxes_.push_back({{0, -1, 0}, {5, 0, 100}});

    BIHBuilder bih(&storage_);
    auto bih_tree = bih(std::move(bboxes_));
    ASSERT_EQ(1, bih_tree.inf_volids.size());
    EXPECT_EQ(LocalVolumeId{0},
              storage_.local_volume_ids[bih_tree.inf_volids[0]]);

    // Test bounding box storage
    auto bbox1 = storage_.bboxes[bih_tree.bboxes[LocalVolumeId{2}]];
    EXPECT_VEC_SOFT_EQ(Real3({1.2, 0., 0.}), bbox1.lower());
    EXPECT_VEC_SOFT_EQ(Real3({2.8, 1., 100.}), bbox1.upper());

    // Test nodes
    auto inner_nodes = bih_tree.inner_nodes;
    auto leaf_nodes = bih_tree.leaf_nodes;
    ASSERT_EQ(3, inner_nodes.size());
    ASSERT_EQ(4, leaf_nodes.size());

    // N0, I0
    {
        auto node = storage_.inner_nodes[inner_nodes[0]];
        ASSERT_FALSE(node.parent);
        EXPECT_EQ(Axis{0}, node.axis);
        EXPECT_EQ(Axis{0}, node.axis);
        EXPECT_SOFT_EQ(2.8, node.bounding_planes[Edge::left].position);
        EXPECT_SOFT_EQ(0, node.bounding_planes[Edge::right].position);
        EXPECT_EQ(1, node.bounding_planes[Edge::left].child.unchecked_get());
        EXPECT_EQ(2, node.bounding_planes[Edge::right].child.unchecked_get());
    }

    // N1, I1
    {
        auto node = storage_.inner_nodes[inner_nodes[1]];
        ASSERT_EQ(BIHNodeId{0}, node.parent);
        EXPECT_EQ(Axis{0}, node.axis);
        EXPECT_EQ(Axis{0}, node.axis);
        EXPECT_SOFT_EQ(1.6, node.bounding_planes[Edge::left].position);
        EXPECT_SOFT_EQ(1.2, node.bounding_planes[Edge::right].position);
        EXPECT_EQ(3, node.bounding_planes[Edge::left].child.unchecked_get());
        EXPECT_EQ(4, node.bounding_planes[Edge::right].child.unchecked_get());
    }

    // N2, I2
    {
        auto node = storage_.inner_nodes[inner_nodes[2]];
        ASSERT_EQ(BIHNodeId{0}, node.parent);
        EXPECT_EQ(Axis{0}, node.axis);
        EXPECT_EQ(Axis{0}, node.axis);
        EXPECT_SOFT_EQ(5, node.bounding_planes[Edge::left].position);
        EXPECT_SOFT_EQ(2.8, node.bounding_planes[Edge::right].position);
        EXPECT_EQ(5, node.bounding_planes[Edge::left].child.unchecked_get());
        EXPECT_EQ(6, node.bounding_planes[Edge::right].child.unchecked_get());
    }

    // N3, L0
    {
        auto node = storage_.leaf_nodes[leaf_nodes[0]];
        ASSERT_EQ(BIHNodeId{1}, node.parent);
        EXPECT_EQ(1, node.vol_ids.size());
        EXPECT_EQ(1,
                  storage_.local_volume_ids[node.vol_ids[0]].unchecked_get());
    }

    // N3, L1
    {
        auto node = storage_.leaf_nodes[leaf_nodes[1]];
        ASSERT_EQ(BIHNodeId{1}, node.parent);
        EXPECT_EQ(1, node.vol_ids.size());
        EXPECT_EQ(2,
                  storage_.local_volume_ids[node.vol_ids[0]].unchecked_get());
    }

    // N5, L2
    {
        auto node = storage_.leaf_nodes[leaf_nodes[2]];
        ASSERT_EQ(BIHNodeId{2}, node.parent);
        EXPECT_EQ(2, node.vol_ids.size());
        EXPECT_EQ(4,
                  storage_.local_volume_ids[node.vol_ids[0]].unchecked_get());
        EXPECT_EQ(5,
                  storage_.local_volume_ids[node.vol_ids[1]].unchecked_get());
    }

    // N6, L3
    {
        auto node = storage_.leaf_nodes[leaf_nodes[3]];
        ASSERT_EQ(BIHNodeId{2}, node.parent);
        EXPECT_EQ(1, node.vol_ids.size());
        EXPECT_EQ(3,
                  storage_.local_volume_ids[node.vol_ids[0]].unchecked_get());
    }
}

//---------------------------------------------------------------------------//
/* Test a 3x4 grid of non-overlapping cuboids.
 *
 *                4 _______________
 *                  | V4 | V8 | V12|
 *                3 |____|____|____|
 *                  | V3 | V7 | V11|
 *            y   2 |____|____|____|
 *                  | V2 | V6 | V10|
 *                1 |____|____|____|
 *                  | V1 | V5 | V9 |
 *                0 |____|____|____|
 *                  0    1    2    3
 *                          x
 *
 * Resultant tree structure in terms of BIHNodeId (N) and volumes (V)
 *
 *                   _______________ N0 ______________
 *                 /                                   \
 *          ___  N1  ___                         ___   N6  ___
 *        /              \                     /                \
 *      N2                N3                 N7                  N8
 *     /   \           /      \             /   \            /       \
 *  N11    N12       N4         N5         N17    N18      N9          N10
 *  V1     V2      /   \      /   \        V3    V4       /  \        /   \
 *                N13  N14   N15   N16                   N19  N20    N21   N22
 *                V5   V6    V9    V10                   V7   V8     V11   V12
 *
 * In terms of BIHInnerNodeIds (I) and BIHLeafNodeIds (L):
 *
 *                   _______________ I0 ______________
 *                 /                                   \
 *          ___  I1  ___                         ___   I6  ___
 *        /              \                     /                \
 *      I2                I3                 I7                 I8
 *     /   \           /      \             /   \            /       \
 *  L0     L1       I4         I5          L6    L7        I9          I10
 *  V1     V2      /   \      /   \        V3    V4       /  \        /   \
 *                L2   L3    L4    L5                    L8   L9     L10   L11
 *                V5   V6    V9    V10                   V7   V8     V11   V12
 *
 *
 *
 * Here, we test only the N1 side for the tree for brevity, as the N6 side is
 * directly analogus.
 */
TEST_F(BIHBuilderTest, grid)
{
    using Edge = BIHInnerNode::Edge;

    bboxes_.push_back({{0, 0, 0}, {1, 1, 100}});
    bboxes_.push_back({{0, 1, 0}, {1, 2, 100}});
    bboxes_.push_back({{0, 2, 0}, {1, 3, 100}});
    bboxes_.push_back({{0, 3, 0}, {1, 4, 100}});

    bboxes_.push_back({{1, 0, 0}, {2, 1, 100}});
    bboxes_.push_back({{1, 1, 0}, {2, 2, 100}});
    bboxes_.push_back({{1, 2, 0}, {2, 3, 100}});
    bboxes_.push_back({{1, 3, 0}, {2, 4, 100}});

    bboxes_.push_back({{2, 0, 0}, {3, 1, 100}});
    bboxes_.push_back({{2, 1, 0}, {3, 2, 100}});
    bboxes_.push_back({{2, 2, 0}, {3, 3, 100}});
    bboxes_.push_back({{2, 3, 0}, {3, 4, 100}});

    BIHBuilder bih(&storage_);
    auto bih_tree = bih(std::move(bboxes_));
    ASSERT_EQ(1, bih_tree.inf_volids.size());
    EXPECT_EQ(LocalVolumeId{0},
              storage_.local_volume_ids[bih_tree.inf_volids[0]]);

    // Test nodes
    auto inner_nodes = bih_tree.inner_nodes;
    auto leaf_nodes = bih_tree.leaf_nodes;
    ASSERT_EQ(11, inner_nodes.size());
    ASSERT_EQ(12, leaf_nodes.size());

    // N0, I0
    {
        auto node = storage_.inner_nodes[inner_nodes[0]];
        ASSERT_FALSE(node.parent);
        EXPECT_EQ(Axis{1}, node.axis);
        EXPECT_EQ(Axis{1}, node.axis);
        EXPECT_SOFT_EQ(2, node.bounding_planes[Edge::left].position);
        EXPECT_SOFT_EQ(2, node.bounding_planes[Edge::right].position);
        EXPECT_EQ(1, node.bounding_planes[Edge::left].child.unchecked_get());
        EXPECT_EQ(6, node.bounding_planes[Edge::right].child.unchecked_get());
    }

    // N1, I1
    {
        auto node = storage_.inner_nodes[inner_nodes[1]];
        ASSERT_EQ(BIHNodeId{0}, node.parent);
        EXPECT_EQ(Axis{0}, node.axis);
        EXPECT_EQ(Axis{0}, node.axis);
        EXPECT_SOFT_EQ(1, node.bounding_planes[Edge::left].position);
        EXPECT_SOFT_EQ(1, node.bounding_planes[Edge::right].position);
        EXPECT_EQ(2, node.bounding_planes[Edge::left].child.unchecked_get());
        EXPECT_EQ(3, node.bounding_planes[Edge::right].child.unchecked_get());
    }

    // N2, I2
    {
        auto node = storage_.inner_nodes[inner_nodes[2]];
        ASSERT_EQ(BIHNodeId{1}, node.parent);
        EXPECT_EQ(Axis{1}, node.axis);
        EXPECT_EQ(Axis{1}, node.axis);
        EXPECT_SOFT_EQ(1, node.bounding_planes[Edge::left].position);
        EXPECT_SOFT_EQ(1, node.bounding_planes[Edge::right].position);
        EXPECT_EQ(11, node.bounding_planes[Edge::left].child.unchecked_get());
        EXPECT_EQ(12, node.bounding_planes[Edge::right].child.unchecked_get());
    }

    // N3, I3
    {
        auto node = storage_.inner_nodes[inner_nodes[3]];
        ASSERT_EQ(BIHNodeId{1}, node.parent);
        EXPECT_EQ(Axis{0}, node.axis);
        EXPECT_EQ(Axis{0}, node.axis);
        EXPECT_SOFT_EQ(2, node.bounding_planes[Edge::left].position);
        EXPECT_SOFT_EQ(2, node.bounding_planes[Edge::right].position);
        EXPECT_EQ(4, node.bounding_planes[Edge::left].child.unchecked_get());
        EXPECT_EQ(5, node.bounding_planes[Edge::right].child.unchecked_get());
    }

    // N4, I4
    {
        auto node = storage_.inner_nodes[inner_nodes[4]];
        ASSERT_EQ(BIHNodeId{3}, node.parent);
        EXPECT_EQ(Axis{1}, node.axis);
        EXPECT_EQ(Axis{1}, node.axis);
        EXPECT_SOFT_EQ(1, node.bounding_planes[Edge::left].position);
        EXPECT_SOFT_EQ(1, node.bounding_planes[Edge::right].position);
        EXPECT_EQ(13, node.bounding_planes[Edge::left].child.unchecked_get());
        EXPECT_EQ(14, node.bounding_planes[Edge::right].child.unchecked_get());
    }

    // N5, I5
    {
        auto node = storage_.inner_nodes[inner_nodes[5]];
        ASSERT_EQ(BIHNodeId{3}, node.parent);
        EXPECT_EQ(Axis{1}, node.axis);
        EXPECT_EQ(Axis{1}, node.axis);
        EXPECT_SOFT_EQ(1, node.bounding_planes[Edge::left].position);
        EXPECT_SOFT_EQ(1, node.bounding_planes[Edge::right].position);
        EXPECT_EQ(15, node.bounding_planes[Edge::left].child.unchecked_get());
        EXPECT_EQ(16, node.bounding_planes[Edge::right].child.unchecked_get());
    }

    // N11, I0
    {
        auto node = storage_.leaf_nodes[leaf_nodes[0]];
        ASSERT_EQ(BIHNodeId{2}, node.parent);
        EXPECT_EQ(1, node.vol_ids.size());
        EXPECT_EQ(1,
                  storage_.local_volume_ids[node.vol_ids[0]].unchecked_get());
    }

    // N12, L1
    {
        auto node = storage_.leaf_nodes[leaf_nodes[1]];
        ASSERT_EQ(BIHNodeId{2}, node.parent);
        EXPECT_EQ(1, node.vol_ids.size());
        EXPECT_EQ(2,
                  storage_.local_volume_ids[node.vol_ids[0]].unchecked_get());
    }

    // N13, L2
    {
        auto node = storage_.leaf_nodes[leaf_nodes[2]];
        ASSERT_EQ(BIHNodeId{4}, node.parent);
        EXPECT_EQ(1, node.vol_ids.size());
        EXPECT_EQ(5,
                  storage_.local_volume_ids[node.vol_ids[0]].unchecked_get());
    }

    // N14, L3
    {
        auto node = storage_.leaf_nodes[leaf_nodes[3]];
        ASSERT_EQ(BIHNodeId{4}, node.parent);
        EXPECT_EQ(1, node.vol_ids.size());
        EXPECT_EQ(6,
                  storage_.local_volume_ids[node.vol_ids[0]].unchecked_get());
    }

    // N15, L4
    {
        auto node = storage_.leaf_nodes[leaf_nodes[4]];
        ASSERT_EQ(BIHNodeId{5}, node.parent);
        EXPECT_EQ(1, node.vol_ids.size());
        EXPECT_EQ(9,
                  storage_.local_volume_ids[node.vol_ids[0]].unchecked_get());
    }

    // N16, L5
    {
        auto node = storage_.leaf_nodes[leaf_nodes[5]];
        ASSERT_EQ(BIHNodeId{5}, node.parent);
        EXPECT_EQ(1, node.vol_ids.size());
        EXPECT_EQ(10,
                  storage_.local_volume_ids[node.vol_ids[0]].unchecked_get());
    }
}

//---------------------------------------------------------------------------//
}  // namespace test
}  // namespace celeritas
