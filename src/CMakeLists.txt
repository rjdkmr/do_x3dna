ENABLE_LANGUAGE(CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

##############################################################
################### CHECKING GROMACS PATH ####################
##############################################################
find_package(GROMACS COMPONENTS "libgromacs")

if (NOT GROMACS_FOUND)
    find_package(GROMACS COMPONENTS "libgromacs_d")
endif(NOT GROMACS_FOUND)

if (NOT GROMACS_FOUND)
    find_package(GROMACS COMPONENTS "libgromacs_mpi")
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

