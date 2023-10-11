# ============== Eclipse CDT Error Handling======
# IF(CMAKE_COMPILER_IS_GNUCC)
#  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fmessage-length=0")
# ENDIF(CMAKE_COMPILER_IS_GNUCC)
# IF(CMAKE_COMPILER_IS_GNUCXX)
#  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmessage-length=0")
# ENDIF(CMAKE_COMPILER_IS_GNUCXX)
# set(CMAKE_VERBOSE_MAKEFILE ON)
# ===============================

##############################################################
################### CHECKING GROMACS PATH ####################
##############################################################
find_package(GROMACS COMPONENTS "libgromacs")
if (NOT GROMACS_FOUND)
	message(STATUS "Continue checking for other Gromacs version library file")
endif(NOT GROMACS_FOUND)

if (NOT GROMACS_FOUND)
	message(FATAL_ERROR "\nGromacs path is not defined.....\nUSE:\n -DGMX_PATH=/path/to/gromacs\n  or  \n export CMAKE_PREFIX_PATH=\"/path/to/gromacs\"\n" )
endif(NOT GROMACS_FOUND)

message(STATUS "Found Gromacs include directory: " ${GROMACS_INCLUDE_DIRS})

if(${GROMACS_VERSION} VERSION_LESS "2023")
	message(FATAL_ERROR, "\n Not implemented for Gromacs version: ${GROMACS_VERSION} \n")
endif()

##############################################################
################### Including GROMACS SRC ####################
##############################################################
if(GMX_SRC)
    if(EXISTS "${GMX_SRC}/src/gromacs")
        set(GROMACS_SRC "${GMX_SRC}/src/")
        message(STATUS "Found GROMACS Source directory for header files: ${GMX_SRC}/src/")
    else()
        message(FATAL_ERROR "\n ${GMX_SRC}/src/gromacs not found in GROMACS SOURCE PATH. \nUSE:\n -DGMX_SRC=/path/to/gromacs/source\n" )
    endif()
else()
    message(FATAL_ERROR "\nGROMACS SOURCE PATH is not defined.....\nUSE:\n -DGMX_SRC=/path/to/gromacs/source\n" )
endif()
##############################################################

##############################################################
# After gromacs2016, it uses c++11, and therefore it must be compiled with c++11 compiler
ENABLE_LANGUAGE(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#include(CheckCXXCompilerFlag)
#CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
#CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
#if(COMPILER_SUPPORTS_CXX11)
#	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
#elseif(COMPILER_SUPPORTS_CXX0X)
#	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
#else()
#	message(FATAL_ERROR "This version of GROMACS requires C++11-compatible standard library. Please use a newer compiler, or a newer standard library, or use the GROMACS 5.1.x release. See the installation guide for details.")
#endif()
##############################################################

FILE(GLOB INCLUDE_DIRS RELATIVE ${GROMACS_SRC} ${GROMACS_SRC}/gromacs/*/include)
FOREACH(INC_DIR ${INCLUDE_DIRS})
  include_directories(${GROMACS_SRC}${INC_DIR})
ENDFOREACH()

include_directories( ${GROMACS_SRC} )
include_directories( ${GROMACS_INCLUDE_DIRS} )
add_definitions( ${GROMACS_DEFINITIONS} )
add_executable(do_x3dna do_x3dna.cpp ExtractData.h ExtractData.cpp)
target_link_libraries(do_x3dna ${GROMACS_LIBRARIES} )
SET_TARGET_PROPERTIES(do_x3dna PROPERTIES LINKER_LANGUAGE CXX)
install(TARGETS do_x3dna DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin)

