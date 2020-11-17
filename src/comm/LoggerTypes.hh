//----------------------------------*-C++-*----------------------------------//
// Copyright 2020 UT-Battelle, LLC, and other Celeritas developers.
// See the top-level COPYRIGHT file for details.
// SPDX-License-Identifier: (Apache-2.0 OR MIT)
//---------------------------------------------------------------------------//
//! \file LoggerTypes.hh
//---------------------------------------------------------------------------//
#pragma once

#include <functional>
#include <string>

namespace celeritas
{
//---------------------------------------------------------------------------//
/*!
 * Enumeration for how important a log message is.
 */
enum class LogLevel
{
    debug,      //!< Debugging messages
    diagnostic, //!< Diagnostics about current program execution
    status,     //!< Program execution status (what stage is beginning)
    info,       //!< Important informational messages
    warning,    //!< Warnings about unusual events
    error,      //!< Something went wrong, but execution continues
    critical    //!< Something went terribly wrong; we're aborting now! Bye!
};

//---------------------------------------------------------------------------//
// Get the plain text equivalent of the log level above
const char* to_cstring(LogLevel);

//---------------------------------------------------------------------------//
//! Stand-in for a more complex class for the "provenance" of data
struct Provenance
{
    std::string file;
    int         line = 0;
};

//! Type for handling a log message
using LogHandler = std::function<void(Provenance, LogLevel, std::string)>;

//---------------------------------------------------------------------------//
} // namespace celeritas
