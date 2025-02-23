#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2009-2011, by the VOTCA Development Team (http://www.votca.org).
# Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
# Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
# and including many others, as listed in the AUTHORS file in the
# top-level source directory and at http://www.gromacs.org.
#
# GROMACS is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.
#
# GROMACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with GROMACS; if not, see
# http://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at http://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out http://www.gromacs.org.

# - Finds parts of GROMACS
# Find the native GROMACS components headers and libraries.
#
#  GROMACS_INCLUDE_DIRS   - where to find GROMACS headers.
#  GROMACS_LIBRARIES      - List of libraries when used by GROMACS.
#  GROMACS_FOUND          - True if all GROMACS components were found.
#  GROMACS_DEFINITIONS    - Extra definies needed by GROMACS
#  GROMACS_PKG            - The name of the pkg-config package needed
#  GROMACS_VERSION        - GROMACS lib interface version
#  GROMACS_MAJOR_VERSION  - GROMACS lib interface major version
#  GROMACS_MINOR_VERSION  - GROMACS lib interface minor version
#  GROMACS_PATCH_LEVEL    - GROMACS lib interface patch level
#  GROMACS_VERSION_STRING - GROMACS lib interface version string (e.g. "4.5.3")
#

########## Look for pkg-config ##########################################
include(GNUInstallDirs)
find_package(PkgConfig)
set(_pkgconfig_path $ENV{PKG_CONFIG_PATH}) # backup original pkg_config_path
#########################################################################

########## To add Path of CMAKE_PREFIX_PATH in PKG_CONFIG_PATH ###########
if (DEFINED CMAKE_PREFIX_PATH_LIST)
  string(REPLACE ":" ";" CMAKE_PREFIX_PATH_LIST $ENV{CMAKE_PREFIX_PATH})
  foreach(_dir ${CMAKE_PREFIX_PATH_LIST})
    if(IS_DIRECTORY "${_dir}/${CMAKE_INSTALL_LIBDIR}/pkgconfig/")
      set(ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:${_dir}/${CMAKE_INSTALL_LIBDIR}/pkgconfig")
    endif()
  endforeach()
endif(DEFINED CMAKE_PREFIX_PATH_LIST)
#########################################################################

####### To add path of gromacs if manually given by user ################
if(DEFINED GMX_PATH)
  if(IS_DIRECTORY "${GMX_PATH}/lib/pkgconfig/")
    set(ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:${GMX_PATH}/lib/pkgconfig")
  endif()
  if(IS_DIRECTORY "${GMX_PATH}/${CMAKE_INSTALL_LIBDIR}/pkgconfig/")
    set(ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:${GMX_PATH}/${CMAKE_INSTALL_LIBDIR}/pkgconfig")
  endif()
endif()
#########################################################################

# extract real library name from component name
if(${GROMACS_FIND_COMPONENTS} MATCHES "^libgromacs(_d|_mpi)?$")
  set(GROMACS_PKG "${GROMACS_FIND_COMPONENTS}")
  string(REGEX REPLACE "^lib(.*)" "\\1" GROMACS_LIBRARY_NAME "${GROMACS_PKG}")
else()
  message(FATAL_ERROR "We do not support finding ${GROMACS_FIND_COMPONENTS}, go and implement it ;-)")
endif()

# Check if static build is requested
if(NOT DEFINED BUILD_STATIC)
  set(BUILD_STATIC OFF)
endif()

if(BUILD_STATIC)
  # Static build is required
  # Unfortunately, GROMACS does not install static libraries by default
  # So we need to find it in the source directory where GROMACS was originally built
  # Assumed that GROMACS was built in ${GMX_SRC}/build.
  # TODO: make this more flexible
  if(EXISTS ${GMX_SRC}/build/lib)  
    find_library(GROMACS_LIBRARY NAMES ${GROMACS_LIBRARY_NAME} HINTS ${GMX_SRC}/build/lib)
    find_library(MUPARSER_LIBRARY NAMES muparser HINTS ${GMX_SRC}/build/lib)
  else()
    message(wARNING "${GMX_SRC}/build/lib does not exist. Static library ${GROMACS_LIBRARY_NAME} not found. Continuing with shared build.")
  endif()

  # IF GROMACS_LIBRARY is found, we need to set some flags
  if(GROMACS_LIBRARY)
    execute_process(COMMAND ${GMX_PATH}/bin/gmx --version COMMAND grep "GROMACS version" 
                    COMMAND awk "{print $3}" OUTPUT_VARIABLE GROMACS_VERSION) # Getting GROMACS version
    find_path(GROMACS_INCLUDE_DIR gromacs HINTS ${GMX_PATH}/include) # assume that include directory is in ${GMX_PATH}/include
    list(APPEND GMX_DEP_LIBRARIES "-fopenmp")
    list(APPEND GMX_DEP_LIBRARIES "${MUPARSER_LIBRARY}")
    list(APPEND GMX_DEP_LIBRARIES "-static")
  else()
    message(wARNING, "Static library ${GROMACS_LIBRARY_NAME} not found in ${GMX_SRC}/build/lib")
  endif()

else(BUILD_STATIC)
  # Static build is not required, so we can use pkg-config to find the shared library
  pkg_check_modules(PC_GROMACS ${GROMACS_PKG})
  set(ENV{PKG_CONFIG_PATH} ${_pkgconfig_path}) # restore original pkg_config_path

  if(PC_GROMACS_FOUND)
    if ("${GROMACS_PKG}" MATCHES "_d$")
      list(APPEND GMX_DEFS "-DGMX_DOUBLE")
    endif("${GROMACS_PKG}" MATCHES "_d$")

    if (PC_GROMACS_CFLAGS_OTHER)
      foreach(DEF ${PC_GROMACS_CFLAGS_OTHER})
        if (${DEF} MATCHES "^-D")
          list(APPEND GMX_DEFS ${DEF})
        endif (${DEF} MATCHES "^-D")
      endforeach(DEF)
      list(REMOVE_DUPLICATES GMX_DEFS)
    endif (PC_GROMACS_CFLAGS_OTHER)
    set(GROMACS_DEFINITIONS "${GMX_DEFS}" CACHE STRING "extra GROMACS definitions")

    find_library(GROMACS_LIBRARY NAMES ${GROMACS_LIBRARY_NAME} HINTS ${PC_GROMACS_LIBRARY_DIRS} ${PC_GROMACS_STATIC_LIBRARY_DIRS} )
    if (GROMACS_LIBRARY)
      if("${GROMACS_LIBRARY}" MATCHES "libgromacs.so")
        if(PC_GROMACS_LIBRARIES)
          list(REMOVE_ITEM PC_GROMACS_LIBRARIES ${GROMACS_LIBRARY_NAME})
          foreach (LIB ${PC_GROMACS_LIBRARIES})
            find_library(GROMACS_${LIB} ${LIB} HINTS ${PC_GROMACS_LIBRARY_DIRS})
            if(GROMACS_${LIB})
              list(APPEND GMX_DEP_LIBRARIES ${GROMACS_${LIB}})
            endif(GROMACS_${LIB})
            unset(GROMACS_${LIB} CACHE)
          endforeach(LIB)
        endif(PC_GROMACS_LIBRARIES)
        if(PC_GROMACS_CFLAGS_OTHER)
          foreach(LIB ${PC_GROMACS_CFLAGS_OTHER})
            if (${LIB} MATCHES "thread")
              find_package(Threads REQUIRED)
              list(APPEND GMX_DEP_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
            endif (${LIB} MATCHES "thread")
          endforeach(LIB)
        endif(PC_GROMACS_CFLAGS_OTHER)

        if( NOT ("${GMX_DEP_LIBRARIES}" MATCHES "pthread") AND "${PC_GROMACS_STATIC_LIBRARIES}" MATCHES "pthread")
          list(APPEND GMX_DEP_LIBRARIES "-lpthread")
        endif(NOT ("${GMX_DEP_LIBRARIES}" MATCHES "pthread") AND "${PC_GROMACS_STATIC_LIBRARIES}" MATCHES "pthread")

        if( NOT ("${GMX_DEP_LIBRARIES}" MATCHES "libm|-lm") )
          list(APPEND GMX_DEP_LIBRARIES "-lm")
        endif(NOT ("${GMX_DEP_LIBRARIES}" MATCHES "libm|-lm") )

        if( NOT ("${GMX_DEP_LIBRARIES}" MATCHES "dl") AND "${PC_GROMACS_STATIC_LIBRARIES}" MATCHES "dl")
          list(APPEND GMX_DEP_LIBRARIES "-ldl")
        endif(NOT ("${GMX_DEP_LIBRARIES}" MATCHES "dl") AND "${PC_GROMACS_STATIC_LIBRARIES}" MATCHES "dl")

        if( NOT ("${GMX_DEP_LIBRARIES}" MATCHES "openmp") AND "${PC_GROMACS_STATIC_LDFLAGS}" MATCHES "openmp")
          list(APPEND GMX_DEP_LIBRARIES "-fopenmp")
        endif(NOT ("${GMX_DEP_LIBRARIES}" MATCHES "openmp") AND "${PC_GROMACS_STATIC_LDFLAGS}" MATCHES "openmp")

        if ( NOT ("${GMX_DEP_LIBRARIES}" MATCHES "fftw3") )
          find_library(FFTW fftw3f HINTS ${FFTW_LIB})
          if (FFTW)
            message(STATUS "Found fftw library: ${FFTW}")
            list(APPEND GMX_DEP_LIBRARIES ${FFTW})
          else(FFTW)
            message(WARNING "\n FFTW library file libfftw3f.so or libfftw3f.a not found... \n In case of FUTURE ERROR relating to fftw,\n USE: -DFFTW_LIB=/path/to/fftw3/lib/\n")
          endif(FFTW)
        endif ( NOT ("${GMX_DEP_LIBRARIES}" MATCHES "fftw3") )

        find_library(ZLIB z HINTS ${ZLIB_PATH})
        if(ZLIB)
          message(STATUS "Found z library: ${ZLIB}")
          list(APPEND GMX_DEP_LIBRARIES ${ZLIB})
        else(ZLIB)
          message(WARNING "libz.a or libz.so not found.\n In case of FUTURE ERROR: undefined reference to `uncompress',\n  Use: -DZLIB_PATH=/path/to/zlib.a(so).\n")
        endif(ZLIB)

      endif("${GROMACS_LIBRARY}" MATCHES "libgromacs.so")
      
      set(GROMACS_VERSION "${PC_GROMACS_VERSION}") # Getting Gromacs version
      find_path(GROMACS_INCLUDE_DIR gromacs HINTS ${PC_GROMACS_INCLUDE_DIRS} ) # Getting Gromacs include directory from pkg-config

      mark_as_advanced(GROMACS_DEFINITIONS)

    endif (GROMACS_LIBRARY)
  endif(PC_GROMACS_FOUND)
endif(BUILD_STATIC)

# Prepare output variables
set(GROMACS_DEP_LIBRARIES ${GMX_DEP_LIBRARIES} CACHE FILEPATH "GROMACS dependency libs (only needed for static (.a) ${GROMACS_LIBRARY}")

# Process GROMACS version
string(REPLACE "." ";" VERSION_LIST ${GROMACS_VERSION})
list(LENGTH VERSION_LIST LEN_GMX_VER_LIST)
list(GET VERSION_LIST 0 GROMACS_MAJOR_VERSION)
if (${LEN_GMX_VER_LIST} MATCHES "2")
  list(GET VERSION_LIST 1 GROMACS_MINOR_VERSION)
endif(${LEN_GMX_VER_LIST} MATCHES "2")
set(GROMACS_VERSION_STRING "${GROMACS_VERSION}")

# Set libraries and the include directories
set(GROMACS_LIBRARIES "${GROMACS_LIBRARY};${GROMACS_DEP_LIBRARIES}" )
set(GROMACS_INCLUDE_DIRS ${GROMACS_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set GROMACS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(GROMACS DEFAULT_MSG GROMACS_LIBRARY GROMACS_INCLUDE_DIR GROMACS_DEP_LIBRARIES)

if(GROMACS_STATIC_LIBRARY)
set(GROMACS_LIBRARY "${GROMACS_STATIC_LIBRARY}")
endif(GROMACS_STATIC_LIBRARY)

mark_as_advanced(GROMACS_INCLUDE_DIR GROMACS_LIBRARY GROMACS_PKG GROMACS_VERSION GROMACS_DEP_LIBRARIES)
