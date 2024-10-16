//----------------------------------*-C++-*----------------------------------//
// Copyright 2022 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file celeritas/global/Stepper.test.cc
//---------------------------------------------------------------------------//
#include "celeritas/global/Stepper.hh"

#include <random>

#include "corecel/Types.hh"
#include "corecel/cont/Range.hh"
#include "celeritas/field/UniformFieldData.hh"
#include "celeritas/global/ActionInterface.hh"
#include "celeritas/global/ActionRegistry.hh"
#include "celeritas/global/alongstep/AlongStepUniformMscAction.hh"
#include "celeritas/phys/PDGNumber.hh"
#include "celeritas/phys/ParticleParams.hh"
#include "celeritas/phys/Primary.hh"
#include "celeritas/random/distribution/IsotropicDistribution.hh"

#include "../SimpleTestBase.hh"
#include "../TestEm15Base.hh"
#include "../TestEm3Base.hh"
#include "DummyAction.hh"
#include "StepperTestBase.hh"
#include "celeritas_test.hh"

using celeritas::units::MevEnergy;

namespace celeritas
{
namespace test
{
//---------------------------------------------------------------------------//
// TEST HARNESS
//---------------------------------------------------------------------------//

#define TestEm3Test TEST_IF_CELERITAS_GEANT(TestEm3Test)
class TestEm3Test : public TestEm3Base, public StepperTestBase
{
  public:
    //! Make 10GeV electrons along +x
    std::vector<Primary> make_primaries(size_type count) const override
    {
        return this->make_primaries_with_energy(count, MevEnergy{10000});
    }

    size_type max_average_steps() const override
    {
        return 100000; // 8 primaries -> ~500k steps, be conservative
    }

    std::vector<Primary>
    make_primaries_with_energy(size_type count, MevEnergy energy) const
    {
        Primary p;
        p.particle_id = this->particle()->find(pdg::electron());
        CELER_ASSERT(p.particle_id);
        p.energy    = energy;
        p.track_id  = TrackId{0};
        p.position  = {-22, 0, 0};
        p.direction = {1, 0, 0};
        p.time      = 0;

        std::vector<Primary> result(count, p);
        for (auto i : range(count))
        {
            result[i].event_id = EventId{i};
        }
        return result;
    }
};

//---------------------------------------------------------------------------//
#define TestEm3MscTest TEST_IF_CELERITAS_GEANT(TestEm3MscTest)
class TestEm3MscTest : public TestEm3Test
{
  public:
    //! Use MSC
    bool enable_msc() const override { return true; }

    //! Make 10MeV electrons along +x
    std::vector<Primary> make_primaries(size_type count) const override
    {
        return this->make_primaries_with_energy(count, MevEnergy{10});
    }

    size_type max_average_steps() const override { return 100; }
};

//---------------------------------------------------------------------------//
#define TestEm3MscNofluctTest TEST_IF_CELERITAS_GEANT(TestEm3MscNofluctTest)
class TestEm3MscNofluctTest : public TestEm3Test
{
  public:
    //! Use MSC
    bool enable_msc() const override { return true; }
    //! Disable energy loss fluctuation
    bool enable_fluctuation() const override { return false; }

    //! Make 10MeV electrons along +x
    std::vector<Primary> make_primaries(size_type count) const override
    {
        return this->make_primaries_with_energy(count, MevEnergy{10});
    }

    size_type max_average_steps() const override { return 100; }
};

//---------------------------------------------------------------------------//
#define TestEm15Test TEST_IF_CELERITAS_GEANT(TestEm15Test)
class TestEm15FieldTest : public TestEm15Base, public StepperTestBase
{
    bool enable_fluctuation() const override { return false; }

    SPConstAction build_along_step() override
    {
        CELER_EXPECT(!this->enable_fluctuation());
        UniformFieldParams field_params;
        field_params.field = {0, 0, 1e-3 * units::tesla};
        auto result        = AlongStepUniformMscAction::from_params(
            *this->physics(), field_params, this->action_reg().get());
        CELER_ENSURE(result);
        CELER_ENSURE(result->has_msc() == this->enable_msc());
        return result;
    }

    //! Make isotropic 10MeV electron/positron mix
    std::vector<Primary> make_primaries(size_type count) const override
    {
        Primary p;
        p.energy   = MevEnergy{10};
        p.position = {0, 0, 0};
        p.time     = 0;
        p.track_id = TrackId{0};

        const Array<ParticleId, 2> particles = {
            this->particle()->find(pdg::electron()),
            this->particle()->find(pdg::positron()),
        };
        CELER_ASSERT(particles[0] && particles[1]);

        std::vector<Primary>    result(count, p);
        IsotropicDistribution<> sample_dir;
        std::mt19937            rng;

        for (auto i : range(count))
        {
            result[i].event_id    = EventId{i};
            result[i].direction   = sample_dir(rng);
            result[i].particle_id = particles[i % particles.size()];
        }
        return result;
    }

    size_type max_average_steps() const override { return 500; }
};

//---------------------------------------------------------------------------//
// TESTEM3
//---------------------------------------------------------------------------//

TEST_F(TestEm3Test, setup)
{
    auto result = this->check_setup();

    static const char* const expected_processes[] = {
        "Compton scattering",
        "Photoelectric effect",
        "Photon annihiliation",
        "Positron annihiliation",
        "Electron/positron ionization",
        "Bremsstrahlung",
    };
    EXPECT_VEC_EQ(expected_processes, result.processes);
    static const char* const expected_actions[] = {
        "pre-step",
        "along-step-general-linear",
        "physics-discrete-select",
        "scat-klein-nishina",
        "photoel-livermore",
        "conv-bethe-heitler",
        "annihil-2-gamma",
        "ioni-moller-bhabha",
        "brems-combined",
        "geo-boundary",
        "dummy-action",
    };
    EXPECT_VEC_EQ(expected_actions, result.actions);
}

TEST_F(TestEm3Test, host)
{
    size_type num_primaries   = 1;
    size_type inits_per_track = 32 * 8;
    size_type num_tracks      = num_primaries * inits_per_track;

    Stepper<MemSpace::host> step(
        this->make_stepper_input(num_tracks, inits_per_track));
    auto result = this->run(step, num_primaries);
    EXPECT_SOFT_NEAR(58000, result.calc_avg_steps_per_primary(), 0.10);

    if (this->is_ci_build())
    {
        EXPECT_EQ(343, result.num_step_iters());
        EXPECT_SOFT_EQ(63490, result.calc_avg_steps_per_primary());
        EXPECT_EQ(255, result.calc_emptying_step());
        EXPECT_EQ(RunResult::StepCount({108, 1416}), result.calc_queue_hwm());
    }
    else if (this->is_wildstyle_build())
    {
        EXPECT_EQ(343, result.num_step_iters());
        EXPECT_SOFT_EQ(63490, result.calc_avg_steps_per_primary());
        EXPECT_EQ(255, result.calc_emptying_step());
        EXPECT_EQ(RunResult::StepCount({108, 1416}), result.calc_queue_hwm());
    }
    else if (this->is_summit_build())
    {
        EXPECT_EQ(323, result.num_step_iters());
        EXPECT_SOFT_EQ(61437, result.calc_avg_steps_per_primary());
        EXPECT_EQ(257, result.calc_emptying_step());
        EXPECT_EQ(RunResult::StepCount({89, 1140}), result.calc_queue_hwm());
    }
    else
    {
        cout << "No output saved for combination of "
             << test::PrintableBuildConf{} << std::endl;
        result.print_expected();

        if (this->strict_testing())
        {
            FAIL() << "Updated stepper results are required for CI tests";
        }
    }

    // Check that callback was called
    EXPECT_EQ(result.active.size(), this->dummy_action().num_execute_host());
    EXPECT_EQ(0, this->dummy_action().num_execute_device());
}

TEST_F(TestEm3Test, TEST_IF_CELER_DEVICE(device))
{
    if (CELERITAS_USE_VECGEOM && this->is_ci_build())
    {
        GTEST_SKIP() << "TODO: TestEm3 + vecgeom crashes on CI";
    }

    size_type num_primaries   = 8;
    size_type inits_per_track = 1024;
    // Num tracks is low enough to hit capacity
    size_type num_tracks = num_primaries * 800;

    Stepper<MemSpace::device> step(
        this->make_stepper_input(num_tracks, inits_per_track));
    auto result = this->run(step, num_primaries);
    EXPECT_SOFT_NEAR(58000, result.calc_avg_steps_per_primary(), 0.10);

    if (this->is_ci_build())
    {
        EXPECT_EQ(218, result.num_step_iters());
        EXPECT_SOFT_EQ(62756.625, result.calc_avg_steps_per_primary());
        EXPECT_EQ(82, result.calc_emptying_step());
        EXPECT_EQ(RunResult::StepCount({75, 1450}), result.calc_queue_hwm());
    }
    else if (this->is_wildstyle_build())
    {
        EXPECT_EQ(218, result.num_step_iters());
        EXPECT_SOFT_EQ(62756.625, result.calc_avg_steps_per_primary());
        EXPECT_EQ(82, result.calc_emptying_step());
        EXPECT_EQ(RunResult::StepCount({75, 1450}), result.calc_queue_hwm());
    }
    else
    {
        cout << "No output saved for combination of "
             << test::PrintableBuildConf{} << std::endl;
        result.print_expected();

        if (this->strict_testing())
        {
            FAIL() << "Updated stepper results are required for CI tests";
        }
    }

    // Check that callback was called
    EXPECT_EQ(result.active.size(), this->dummy_action().num_execute_device());
    EXPECT_EQ(0, this->dummy_action().num_execute_host());
}

//---------------------------------------------------------------------------//
// TESTEM3_MSC
//---------------------------------------------------------------------------//

TEST_F(TestEm3MscTest, setup)
{
    auto result = this->check_setup();

    static const char* const expected_processes[] = {
        "Compton scattering",
        "Photoelectric effect",
        "Photon annihiliation",
        "Positron annihiliation",
        "Electron/positron ionization",
        "Bremsstrahlung",
        "Multiple scattering",
    };
    EXPECT_VEC_EQ(expected_processes, result.processes);
    static const char* const expected_actions[] = {
        "pre-step",
        "along-step-general-linear",
        "physics-discrete-select",
        "scat-klein-nishina",
        "photoel-livermore",
        "conv-bethe-heitler",
        "annihil-2-gamma",
        "ioni-moller-bhabha",
        "brems-combined",
        "msc-urban",
        "geo-boundary",
        "dummy-action",
    };
    EXPECT_VEC_EQ(expected_actions, result.actions);
}

TEST_F(TestEm3MscTest, host)
{
    size_type num_primaries   = 8;
    size_type inits_per_track = 32 * 8;
    size_type num_tracks      = num_primaries * inits_per_track;

    Stepper<MemSpace::host> step(
        this->make_stepper_input(num_tracks, inits_per_track));
    auto result = this->run(step, num_primaries);
    EXPECT_SOFT_NEAR(30.5, result.calc_avg_steps_per_primary(), 0.25);

    if (this->is_ci_build() || this->is_wildstyle_build())
    {
        EXPECT_EQ(30, result.num_step_iters());
        EXPECT_SOFT_EQ(30.625, result.calc_avg_steps_per_primary());
        EXPECT_EQ(10, result.calc_emptying_step());
        EXPECT_EQ(RunResult::StepCount({8, 6}), result.calc_queue_hwm());
    }
    else
    {
        cout << "No output saved for combination of "
             << test::PrintableBuildConf{} << std::endl;
        result.print_expected();

        if (this->strict_testing())
        {
            FAIL() << "Updated stepper results are required for CI tests";
        }
    }
}

TEST_F(TestEm3MscTest, TEST_IF_CELER_DEVICE(device))
{
    size_type num_primaries   = 8;
    size_type inits_per_track = 512;
    size_type num_tracks      = 1024;

    Stepper<MemSpace::device> step(
        this->make_stepper_input(num_tracks, inits_per_track));
    auto result = this->run(step, num_primaries);

    if (this->is_ci_build())
    {
        if (CELERITAS_USE_VECGEOM)
        {
            EXPECT_EQ(64, result.num_step_iters());
            EXPECT_SOFT_EQ(62.5, result.calc_avg_steps_per_primary());
        }
        else
        {
            EXPECT_EQ(63, result.num_step_iters());
            EXPECT_SOFT_EQ(62.375, result.calc_avg_steps_per_primary());
        }
        EXPECT_EQ(8, result.calc_emptying_step());
        EXPECT_EQ(RunResult::StepCount({6, 7}), result.calc_queue_hwm());
    }
    else
    {
        cout << "No output saved for combination of "
             << test::PrintableBuildConf{} << std::endl;
        result.print_expected();

        if (this->strict_testing())
        {
            FAIL() << "Updated stepper results are required for CI tests";
        }
    }
}

//---------------------------------------------------------------------------//
// TESTEM3_MSC_NOFLUCT
//---------------------------------------------------------------------------//

TEST_F(TestEm3MscNofluctTest, host)
{
    size_type num_primaries   = 8;
    size_type inits_per_track = 32 * 8;
    size_type num_tracks      = num_primaries * inits_per_track;

    Stepper<MemSpace::host> step(
        this->make_stepper_input(num_tracks, inits_per_track));
    auto result = this->run(step, num_primaries);
    EXPECT_SOFT_NEAR(55, result.calc_avg_steps_per_primary(), 0.50);

    if (this->is_ci_build())
    {
        EXPECT_EQ(71, result.num_step_iters());
        EXPECT_SOFT_EQ(57.125, result.calc_avg_steps_per_primary());
        EXPECT_EQ(8, result.calc_emptying_step());
        EXPECT_EQ(RunResult::StepCount({4, 5}), result.calc_queue_hwm());
    }
    else
    {
        cout << "No output saved for combination of "
             << test::PrintableBuildConf{} << std::endl;
        result.print_expected();

        if (this->strict_testing())
        {
            FAIL() << "Updated stepper results are required for CI tests";
        }
    }
}

TEST_F(TestEm3MscNofluctTest, TEST_IF_CELER_DEVICE(device))
{
    if (CELERITAS_USE_VECGEOM && this->is_ci_build())
    {
        GTEST_SKIP() << "TODO: TestEm3 + vecgeom crashes on CI";
    }

    size_type num_primaries   = 8;
    size_type inits_per_track = 512;
    size_type num_tracks      = 1024;

    Stepper<MemSpace::device> step(
        this->make_stepper_input(num_tracks, inits_per_track));
    auto result = this->run(step, num_primaries);

    if (this->is_ci_build())
    {
        EXPECT_EQ(38, result.num_step_iters());
        EXPECT_SOFT_EQ(44.75, result.calc_avg_steps_per_primary());
        EXPECT_EQ(11, result.calc_emptying_step());
        EXPECT_EQ(RunResult::StepCount({10, 5}), result.calc_queue_hwm());
    }
    else
    {
        cout << "No output saved for combination of "
             << test::PrintableBuildConf{} << std::endl;
        result.print_expected();

        if (this->strict_testing())
        {
            FAIL() << "Updated stepper results are required for CI tests";
        }
    }
}

//---------------------------------------------------------------------------//
// TESTEM15_MSC_FIELD
//---------------------------------------------------------------------------//

TEST_F(TestEm15FieldTest, setup)
{
    auto result = this->check_setup();

    static const char* const expected_processes[] = {
        "Compton scattering",
        "Photoelectric effect",
        "Photon annihiliation",
        "Positron annihiliation",
        "Electron/positron ionization",
        "Bremsstrahlung",
        "Multiple scattering",
    };
    EXPECT_VEC_EQ(expected_processes, result.processes);
    static const char* const expected_actions[] = {
        "pre-step",
        "along-step-uniform-msc",
        "physics-discrete-select",
        "scat-klein-nishina",
        "photoel-livermore",
        "conv-bethe-heitler",
        "annihil-2-gamma",
        "ioni-moller-bhabha",
        "brems-sb",
        "brems-rel",
        "msc-urban",
        "geo-boundary",
        "dummy-action",
    };
    EXPECT_VEC_EQ(expected_actions, result.actions);
}

TEST_F(TestEm15FieldTest, host)
{
    size_type num_primaries   = 4;
    size_type inits_per_track = 32 * 8;
    size_type num_tracks      = num_primaries * inits_per_track;

    Stepper<MemSpace::host> step(
        this->make_stepper_input(num_tracks, inits_per_track));
    auto result = this->run(step, num_primaries);
    EXPECT_SOFT_NEAR(35, result.calc_avg_steps_per_primary(), 0.50);

    if (this->is_ci_build() || this->is_summit_build()
        || this->is_wildstyle_build())
    {
        EXPECT_EQ(14, result.num_step_iters());
        EXPECT_SOFT_EQ(35, result.calc_avg_steps_per_primary());
        EXPECT_EQ(6, result.calc_emptying_step());
        EXPECT_EQ(RunResult::StepCount({4, 7}), result.calc_queue_hwm());
    }
    else
    {
        cout << "No output saved for combination of "
             << test::PrintableBuildConf{} << std::endl;
        result.print_expected();

        if (this->strict_testing())
        {
            FAIL() << "Updated stepper results are required for CI tests";
        }
    }
}

TEST_F(TestEm15FieldTest, TEST_IF_CELER_DEVICE(device))
{
    size_type num_primaries   = 8;
    size_type inits_per_track = 512;
    size_type num_tracks      = 1024;

    Stepper<MemSpace::device> step(
        this->make_stepper_input(num_tracks, inits_per_track));
    auto result = this->run(step, num_primaries);

    if (this->is_ci_build() || this->is_summit_build()
        || this->is_wildstyle_build())
    {
        EXPECT_EQ(14, result.num_step_iters());
        EXPECT_SOFT_EQ(29.75, result.calc_avg_steps_per_primary());
        EXPECT_EQ(5, result.calc_emptying_step());
        EXPECT_EQ(RunResult::StepCount({2, 11}), result.calc_queue_hwm());
    }
    else
    {
        cout << "No output saved for combination of "
             << test::PrintableBuildConf{} << std::endl;
        result.print_expected();

        if (this->strict_testing())
        {
            FAIL() << "Updated stepper results are required for CI tests";
        }
    }
}

//---------------------------------------------------------------------------//
} // namespace test
} // namespace celeritas
