#----------------------------------*-CMake-*----------------------------------#
# Copyright 2020-2022 UT-Battelle, LLC, and other Celeritas developers.
# See the top-level COPYRIGHT file for details.
# SPDX-License-Identifier: (Apache-2.0 OR MIT)
#[=======================================================================[.rst:

CeleritasUtils
--------------

CMake configuration utility functions for Celeritas.

.. command:: celeritas_optional_language

  Add an configurable cache option ``CELERITAS_USE_<lang>`` that defaults to
  checking whether the language is available.

    celeritas_optional_language(<lang>)

.. command:: celeritas_optional_package

  Add an configurable cache option ``CELERITAS_USE_<package>`` that searches for
  the package to decide its default value.

    celeritas_optional_package(<package> [<find_package>] <docstring>)

  This won't be used for all Celeritas options or even all external dependent
  packages. If given, the ``<find_package>`` package name will searched for
  instead of ``<package>``.

.. command:: celeritas_check_python_module

   Determine whether a given Python module is available with the current
   environment. ::

     celeritas_check_python_module(<variable> <module>)

   ``<variable>``
     Variable name that will be set to whether the module exists

   ``<module>``
     Python module name, e.g. "numpy" or "scipy.linalg"

   Note that because this function caches the Python script result to save
   reconfigure time (or when multiple scripts check for the same module),
   changing the Python executable or installed modules may mean
   having to delete or modify your CMakeCache.txt file.

   Example::

      celeritas_check_python_module(has_numpy "numpy")

#]=======================================================================]
include_guard(GLOBAL)

include(CheckLanguage)

#-----------------------------------------------------------------------------#

function(celeritas_optional_language lang)
  set(_var "CELERITAS_USE_${lang}")
  if(DEFINED "${_var}")
    set(_val "${_var}")
  else()
    check_language(${lang})
    set(_val OFF)
    if(CMAKE_${lang}_COMPILER)
      set(_val ON)
    endif()
    message(STATUS "Set ${_var}=${_val} based on compiler availability")
  endif()

  option("${_var}" "Enable the ${lang} language" "${_val}" )
endfunction()

#-----------------------------------------------------------------------------#

function(celeritas_to_onoff varname)
  if(ARGC GREATER 1 AND ARGV1)
    set(${varname} ON PARENT_SCOPE)
  else()
    set(${varname} OFF PARENT_SCOPE)
  endif()
endfunction()

#-----------------------------------------------------------------------------#

# Note: this is a macro so that `find_package` variables stay in the global
# scope.
macro(celeritas_optional_package package)
  if("${ARGC}" EQUAL 2)
    set(_findpkg "${package}")
    set(_docstring "${ARGV1}")
  else()
    set(_findpkg "${ARGV1}")
    set(_docstring "${ARGV2}")
  endif()

  set(_var "CELERITAS_USE_${package}")
  if(DEFINED "${_var}")
    set(_val "${_var}")
  else()
    find_package(${_findpkg} QUIET)
    list(GET _findpkg 0 _findpkg)
    celeritas_to_onoff(_val ${${_findpkg}_FOUND})
    message(STATUS "Set ${_var}=${_val} based on package availability")
  endif()

  option("${_var}" "${_docstring}" "${_val}")
endmacro()

#-----------------------------------------------------------------------------#

function(celeritas_check_python_module varname module)
  set(_cache_name CELERITAS_CHECK_PYTHON_MODULE_${module})
  if(DEFINED ${_cache_name})
    # We've already checked for this module
    set(_found "${${_cache_name}}")
  else()
    message(STATUS "Check Python module ${module}")
    set(_cmd
      "${CMAKE_COMMAND}" -E env "PYTHONPATH=${CELERITAS_PYTHONPATH}"
      "${Python_EXECUTABLE}" -c "import ${module}"
    )
    execute_process(COMMAND
      ${_cmd}
      RESULT_VARIABLE _result
      ERROR_QUIET # hide error message if module unavailable
    )
    # Note: use JSON-compatible T/F representation
    if(_result)
      set(_msg "not found")
      set(_found false)
    else()
      set(_msg "found")
      set(_found true)
    endif()
    message(STATUS "Check Python module ${module} -- ${_msg}")
    set(${_cache_name} "${_found}" CACHE INTERNAL
      "Whether Python module ${module} is available")
  endif()

  # Save outgoing variable
  set(${varname} "${_found}" PARENT_SCOPE)
endfunction()

#-----------------------------------------------------------------------------#
