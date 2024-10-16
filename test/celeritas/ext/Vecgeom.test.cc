//----------------------------------*-C++-*----------------------------------//
// Copyright 2020-2022 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celeritas/ext/Vecgeom.test.cc
//---------------------------------------------------------------------------//
#include "Vecgeom.test.hh"

#include "corecel/cont/ArrayIO.hh"
#include "corecel/data/CollectionStateStore.hh"
#include "corecel/io/Repr.hh"
#include "corecel/math/NumericLimits.hh"
#include "corecel/sys/Device.hh"
#include "celeritas/GlobalGeoTestBase.hh"
#include "celeritas/GlobalTestBase.hh"
#include "celeritas/ext/LoadGdml.hh"
#include "celeritas/ext/VecgeomData.hh"
#include "celeritas/ext/VecgeomParams.hh"
#include "celeritas/ext/VecgeomTrackView.hh"

#include "celeritas_test.hh"

namespace celeritas
{
namespace test
{
//---------------------------------------------------------------------------//

// Since VecGeom is currently CUDA-only, we cannot use the TEST_IF_CELER_DEVICE
// macro (which also allows HIP).
#if CELERITAS_USE_CUDA
#    define TEST_IF_CELERITAS_CUDA(name) name
#else
#    define TEST_IF_CELERITAS_CUDA(name) DISABLED_##name
#endif

// Always compile and sometimes disable tests that require Geant4
#if CELERITAS_USE_GEANT4
#    define TEST_IF_CELERITAS_GEANT(name) name
#else
#    define TEST_IF_CELERITAS_GEANT(name) DISABLED_##name
#endif

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
class VecgeomTestBase : virtual public GlobalTestBase
{
  public:
    //!@{
    using HostStateStore
        = CollectionStateStore<VecgeomStateData, MemSpace::host>;
    //!@}

    struct TrackingResult
    {
        std::vector<std::string> volumes;
        std::vector<real_type>   distances;

        void print_expected();
    };

  public:
    //! Construct host state (and load geometry) during steup
    void SetUp() override
    {
        host_state = HostStateStore(this->geometry()->host_ref(), 1);
    }

    //! Create a host track view
    VecgeomTrackView make_geo_track_view()
    {
        return VecgeomTrackView(
            this->geometry()->host_ref(), host_state.ref(), ThreadId(0));
    }

    //! Find linear segments until outside
    TrackingResult track(const Real3& pos, const Real3& dir);

  protected:
    SPConstParticle    build_particle() final { CELER_ASSERT_UNREACHABLE(); }
    SPConstCutoff      build_cutoff() final { CELER_ASSERT_UNREACHABLE(); }
    SPConstPhysics     build_physics() final { CELER_ASSERT_UNREACHABLE(); }
    SPConstAction      build_along_step() final { CELER_ASSERT_UNREACHABLE(); }
    SPConstMaterial    build_material() final { CELER_ASSERT_UNREACHABLE(); }
    SPConstGeoMaterial build_geomaterial() final
    {
        CELER_ASSERT_UNREACHABLE();
    }

  private:
    HostStateStore host_state;
};

auto VecgeomTestBase::track(const Real3& pos, const Real3& dir)
    -> TrackingResult
{
    const auto& params = *this->geometry();

    TrackingResult result;

    VecgeomTrackView geo = this->make_geo_track_view();
    geo                  = {pos, dir};

    if (geo.is_outside())
    {
        // Initial step is outside but may approach insidfe
        result.volumes.push_back("[OUTSIDE]");
        auto next = geo.find_next_step();
        result.distances.push_back(next.distance);
        if (next.boundary)
        {
            geo.move_to_boundary();
            geo.cross_boundary();
            EXPECT_TRUE(geo.is_on_boundary());
        }
    }

    while (!geo.is_outside())
    {
        result.volumes.push_back(params.id_to_label(geo.volume_id()).name);
        auto next = geo.find_next_step();
        result.distances.push_back(next.distance);
        if (!next.boundary)
        {
            // Failure to find the next boundary while inside the geometry
            ADD_FAILURE();
            result.volumes.push_back("[NO INTERCEPT]");
            break;
        }
        geo.move_to_boundary();
        geo.cross_boundary();
    }

    return result;
}

void VecgeomTestBase::TrackingResult::print_expected()
{
    cout << "/*** ADD THE FOLLOWING UNIT TEST CODE ***/\n"
         << "static const char* const expected_volumes[] = "
         << repr(this->volumes) << ";\n"
         << "EXPECT_VEC_EQ(expected_volumes, result.volumes);\n"
         << "static const real_type expected_distances[] = "
         << repr(this->distances) << ";\n"
         << "EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);\n"
         << "/*** END CODE ***/\n";
}

//---------------------------------------------------------------------------//

class FourLevelsTest : public VecgeomTestBase, public GlobalGeoTestBase
{
  public:
    const char* geometry_basename() const final { return "four-levels"; }
};

//---------------------------------------------------------------------------//

TEST_F(FourLevelsTest, accessors)
{
    const auto& geom = *this->geometry();
    EXPECT_EQ(4, geom.num_volumes());
    EXPECT_EQ(4, geom.max_depth());

    EXPECT_EQ("Shape2", geom.id_to_label(VolumeId{0}).name);
    EXPECT_EQ("Shape1", geom.id_to_label(VolumeId{1}).name);
    EXPECT_EQ("Envelope", geom.id_to_label(VolumeId{2}).name);
    EXPECT_EQ("World", geom.id_to_label(VolumeId{3}).name);
    EXPECT_EQ(Label("World", "0xdeadbeef"), geom.id_to_label(VolumeId{3}));
}

//---------------------------------------------------------------------------//

TEST_F(FourLevelsTest, detailed_track)
{
    VecgeomTrackView geo = this->make_geo_track_view();
    geo                  = GeoTrackInitializer{{-10, -10, -10}, {1, 0, 0}};
    EXPECT_EQ(VolumeId{0}, geo.volume_id());
    EXPECT_FALSE(geo.is_on_boundary());

    // Check for surfaces up to a distance of 4 units away
    auto next = geo.find_next_step(4.0);
    EXPECT_SOFT_EQ(4.0, next.distance);
    EXPECT_FALSE(next.boundary);
    next = geo.find_next_step(4.0);
    EXPECT_SOFT_EQ(4.0, next.distance);
    EXPECT_FALSE(next.boundary);
    geo.move_internal(3.5);
    EXPECT_FALSE(geo.is_on_boundary());

    // Find one a bit further, then cross it
    next = geo.find_next_step(4.0);
    EXPECT_SOFT_EQ(1.5, next.distance);
    EXPECT_TRUE(next.boundary);
    geo.move_to_boundary();
    EXPECT_EQ(VolumeId{0}, geo.volume_id());
    geo.cross_boundary();
    EXPECT_EQ(VolumeId{1}, geo.volume_id());
    EXPECT_TRUE(geo.is_on_boundary());

    // Find the next boundary and make sure that nearer distances aren't
    // accepted
    next = geo.find_next_step();
    EXPECT_SOFT_EQ(1.0, next.distance);
    EXPECT_TRUE(next.boundary);
    EXPECT_TRUE(geo.is_on_boundary());
    next = geo.find_next_step(0.5);
    EXPECT_SOFT_EQ(0.5, next.distance);
    EXPECT_FALSE(next.boundary);

    {
        SCOPED_TRACE("outside in");
        geo = GeoTrackInitializer{{-25, 6.5, 6.5}, {1, 0, 0}};
        EXPECT_TRUE(geo.is_outside());

        next = geo.find_next_step(0.5);
        EXPECT_SOFT_EQ(0.5, next.distance);
        EXPECT_FALSE(next.boundary);

        next = geo.find_next_step(2);
        EXPECT_SOFT_EQ(1.0, next.distance);
        EXPECT_TRUE(next.boundary);

        geo.move_to_boundary();
        EXPECT_TRUE(geo.is_outside());
        geo.cross_boundary();
        EXPECT_FALSE(geo.is_outside());
        EXPECT_EQ(VolumeId{3}, geo.volume_id());
    }
    {
        SCOPED_TRACE("inside out");
        geo = GeoTrackInitializer{{-23.5, 6.5, 6.5}, {-1, 0, 0}};
        EXPECT_FALSE(geo.is_outside());
        EXPECT_EQ(VolumeId{3}, geo.volume_id());

        next = geo.find_next_step(2);
        EXPECT_SOFT_EQ(0.5, next.distance);
        EXPECT_TRUE(next.boundary);

        geo.move_to_boundary();
        EXPECT_FALSE(geo.is_outside());
        geo.cross_boundary();
        EXPECT_TRUE(geo.is_outside());

        next = geo.find_next_step(2);
        EXPECT_SOFT_EQ(2, next.distance);
        EXPECT_FALSE(next.boundary);

        next = geo.find_next_step();
        EXPECT_GT(next.distance, 1e10);
        EXPECT_FALSE(next.boundary);
    }
}

//---------------------------------------------------------------------------//

TEST_F(FourLevelsTest, tracking)
{
    {
        SCOPED_TRACE("Rightward");
        auto result = this->track({-10, -10, -10}, {1, 0, 0});
        // result.print_expected();
        static const char* const expected_volumes[] = {"Shape2",
                                                       "Shape1",
                                                       "Envelope",
                                                       "World",
                                                       "Envelope",
                                                       "Shape1",
                                                       "Shape2",
                                                       "Shape1",
                                                       "Envelope",
                                                       "World"};
        EXPECT_VEC_EQ(expected_volumes, result.volumes);
        static const real_type expected_distances[]
            = {5, 1, 1, 6, 1, 1, 10, 1, 1, 7};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }
    {
        SCOPED_TRACE("From outside edge");
        auto result = this->track({-24, 10., 10.}, {1, 0, 0});
        static const char* const expected_volumes[] = {"[OUTSIDE]",
                                                       "World",
                                                       "Envelope",
                                                       "Shape1",
                                                       "Shape2",
                                                       "Shape1",
                                                       "Envelope",
                                                       "World",
                                                       "Envelope",
                                                       "Shape1",
                                                       "Shape2",
                                                       "Shape1",
                                                       "Envelope",
                                                       "World"};
        EXPECT_VEC_EQ(expected_volumes, result.volumes);
        static const real_type expected_distances[]
            = {1e-13, 7.0 - 1e-13, 1, 1, 10, 1, 1, 6, 1, 1, 10, 1, 1, 7};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }
    {
        SCOPED_TRACE("Leaving world");
        auto result = this->track({-10, 10, 10}, {0, 1, 0});
        static const char* const expected_volumes[]
            = {"Shape2", "Shape1", "Envelope", "World"};
        EXPECT_VEC_EQ(expected_volumes, result.volumes);
        static const real_type expected_distances[] = {5, 1, 2, 6};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }
    {
        SCOPED_TRACE("Upward");
        auto result = this->track({-10, 10, 10}, {0, 0, 1});
        static const char* const expected_volumes[]
            = {"Shape2", "Shape1", "Envelope", "World"};
        EXPECT_VEC_EQ(expected_volumes, result.volumes);
        static const real_type expected_distances[] = {5, 1, 3, 5};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }
    {
        // Formerly in linear propagator test, used to fail
        SCOPED_TRACE("From just outside world");
        auto result = this->track({-24, 6.5, 6.5}, {1, 0, 0});
        static const char* const expected_volumes[] = {"[OUTSIDE]",
                                                       "World",
                                                       "Envelope",
                                                       "Shape1",
                                                       "Shape2",
                                                       "Shape1",
                                                       "Envelope",
                                                       "World",
                                                       "Envelope",
                                                       "Shape1",
                                                       "Shape2",
                                                       "Shape1",
                                                       "Envelope",
                                                       "World"};
        EXPECT_VEC_EQ(expected_volumes, result.volumes);
        static const real_type expected_distances[] = {1e-13,
                                                       6.9999999999999,
                                                       1,
                                                       5.2928932188135,
                                                       1.4142135623731,
                                                       5.2928932188135,
                                                       1,
                                                       6,
                                                       1,
                                                       5.2928932188135,
                                                       1.4142135623731,
                                                       5.2928932188135,
                                                       1,
                                                       7};
        EXPECT_VEC_SOFT_EQ(expected_distances, result.distances);
    }
}

//---------------------------------------------------------------------------//

TEST_F(FourLevelsTest, safety)
{
    VecgeomTrackView       geo = this->make_geo_track_view();
    std::vector<real_type> safeties;

    for (auto i : range(11))
    {
        real_type r = 2.0 * i;
        geo         = {{r, r, r}, {1, 0, 0}};

        if (!geo.is_outside())
        {
            safeties.push_back(geo.find_safety());
        }
    }

    static const real_type expected_safeties[] = {3,
                                                  1,
                                                  0,
                                                  1.92820323027551,
                                                  1.53589838486225,
                                                  5,
                                                  1.53589838486225,
                                                  1.92820323027551,
                                                  0,
                                                  1,
                                                  3};
    EXPECT_VEC_SOFT_EQ(expected_safeties, safeties);
}

//---------------------------------------------------------------------------//

TEST_F(FourLevelsTest, TEST_IF_CELERITAS_CUDA(device))
{
    using StateStore = CollectionStateStore<VecgeomStateData, MemSpace::device>;

    // Set up test input
    VGGTestInput input;
    input.init = {{{10, 10, 10}, {1, 0, 0}},
                  {{10, 10, -10}, {1, 0, 0}},
                  {{10, -10, 10}, {1, 0, 0}},
                  {{10, -10, -10}, {1, 0, 0}},
                  {{-10, 10, 10}, {-1, 0, 0}},
                  {{-10, 10, -10}, {-1, 0, 0}},
                  {{-10, -10, 10}, {-1, 0, 0}},
                  {{-10, -10, -10}, {-1, 0, 0}}};
    StateStore device_states(this->geometry()->host_ref(), input.init.size());
    input.max_segments = 5;
    input.params       = this->geometry()->device_ref();
    input.state        = device_states.ref();

    // Run kernel
    auto output = vgg_test(input);

    static const int expected_ids[]
        = {1, 2, 3, -2, -3, 1, 2, 3, -2, -3, 1, 2, 3, -2, -3, 1, 2, 3, -2, -3,
           1, 2, 3, -2, -3, 1, 2, 3, -2, -3, 1, 2, 3, -2, -3, 1, 2, 3, -2, -3};

    static const double expected_distances[]
        = {5, 1, 1, 7, -3, 5, 1, 1, 7, -3, 5, 1, 1, 7, -3, 5, 1, 1, 7, -3,
           5, 1, 1, 7, -3, 5, 1, 1, 7, -3, 5, 1, 1, 7, -3, 5, 1, 1, 7, -3};

    // Check results
    EXPECT_VEC_EQ(expected_ids, output.ids);
    EXPECT_VEC_SOFT_EQ(expected_distances, output.distances);
}

//---------------------------------------------------------------------------//
// CONSTRUCT FROM GEANT4 (TODO)
//---------------------------------------------------------------------------//

#define GeantBuilderTest TEST_IF_CELERITAS_GEANT(GeantBuilderTest)
class GeantBuilderTest : public VecgeomTestBase, virtual public GlobalTestBase
{
  public:
    static void SetUpTestCase()
    {
        // Make sure existing VecGeom geometry has been cleared
        test::GlobalGeoTestBase::reset_geometry();
    }

    void SetUp() override
    {
        VecgeomTestBase::SetUp();
        world_volume_
            = load_gdml(this->test_data_path("celeritas", "four-levels.gdml"));
    }

    SPConstGeo build_geometry() override
    {
        CELER_NOT_IMPLEMENTED("build_geometry");
        // return std::make_shared<VecgeomParams>(world_volume.get());
    }

  private:
    UPG4PhysicalVolume world_volume_;
};
//---------------------------------------------------------------------------//
} // namespace test
} // namespace celeritas
