/*----------------------------------*-C-*------------------------------------*
 * Copyright 2020-2022 UT-Battelle, LLC, and other Celeritas developers.
 * See the top-level COPYRIGHT file for details.
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *---------------------------------------------------------------------------*/
/*! \file celeritas_version.h
 * Version metadata for Celeritas.
 *---------------------------------------------------------------------------*/
#ifndef celeritas_version_h
#define celeritas_version_h

#ifdef __cplusplus
extern "C" {
#endif

extern const char celeritas_version[];     /*!< Full version string */
extern const int  celeritas_version_major; /*!< First version component */
extern const int  celeritas_version_minor; /*!< Second version component */
extern const int  celeritas_version_patch; /*!< Third version component */

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* celeritas_version_h */
