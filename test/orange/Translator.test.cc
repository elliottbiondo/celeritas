//----------------------------------*-C++-*----------------------------------//
// Copyright 2021-2022 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file orange/Translator.test.cc
//---------------------------------------------------------------------------//
#include "orange/Translator.hh"

#include "celeritas_test.hh"

namespace celeritas
{
namespace test
{
//---------------------------------------------------------------------------//
using TranslatorTest = Test;

TEST_F(TranslatorTest, basic)
{
    Translator trans;
    Real3      point{0.1, 0.2, 0.3};

    trans(point, {1, 2, 3});

    EXPECT_VEC_SOFT_EQ((Real3{1.1, 2.2, 3.3}), point);
}

//---------------------------------------------------------------------------//
} // namespace test
} // namespace celeritas
