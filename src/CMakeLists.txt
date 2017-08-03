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
	find_package(GROMACS COMPONENTS "libgmx")
endif(NOT GROMACS_FOUND)

if (NOT GROMACS_FOUND)
	message(FATAL_ERROR "\nGromacs path is not defined.....\nUSE:\n -DGMX_PATH=/path/to/gromacs\n  or  \n export CMAKE_PREFIX_PATH=\"/path/to/gromacs\"\n" )
endif(NOT GROMACS_FOUND)

message(STATUS "Found Gromacs include directory: " ${GROMACS_INCLUDE_DIRS})

if(${GROMACS_VERSION} VERSION_LESS "4.5.0")
	message(FATAL_ERROR, "\n Not implemented for Gromacs version: ${GROMACS_VERSION} \n")
endif()

if(${GROMACS_VERSION} VERSION_GREATER "4.5.0" AND ${GROMACS_VERSION} VERSION_LESS "5.0")
	set(HAVE_GROMACS45 TRUE)
	ADD_DEFINITIONS("-DHAVE_GROMACS45")
endif()

if(${GROMACS_VERSION} VERSION_GREATER "5.0" AND ${GROMACS_VERSION} VERSION_LESS "5.1")
	set(HAVE_GROMACS50 TRUE)
	ADD_DEFINITIONS("-DHAVE_GROMACS50")
  ENABLE_LANGUAGE(CXX)
endif(${GROMACS_VERSION} VERSION_GREATER "5.0" AND ${GROMACS_VERSION} VERSION_LESS "5.1")

if(${GROMACS_VERSION} VERSION_GREATER "5.1" AND ${GROMACS_VERSION} VERSION_LESS "5.2" )
	set(HAVE_GROMACS51 TRUE)
	ADD_DEFINITIONS("-DHAVE_GROMACS51")
  ENABLE_LANGUAGE(CXX)
endif(${GROMACS_VERSION} VERSION_GREATER "5.1" AND ${GROMACS_VERSION} VERSION_LESS "5.2")

if(${GROMACS_VERSION} VERSION_GREATER "2016" AND ${GROMACS_VERSION} VERSION_LESS "2017" )
	set(HAVE_GROMACS2016 TRUE)
	ADD_DEFINITIONS("-DHAVE_GROMACS2016")
  ENABLE_LANGUAGE(CXX)
endif(${GROMACS_VERSION} VERSION_GREATER "2016" AND ${GROMACS_VERSION} VERSION_LESS "2017")

include_directories( ${GROMACS_INCLUDE_DIRS} )
add_definitions( ${GROMACS_DEFINITIONS} )
##############################################################

##############################################################
# After gromacs2016, it uses c++11, and therefore it must be compiled with c++11 compiler
if( HAVE_GROMACS2016 )
	include(CheckCXXCompilerFlag)
	CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
	CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
	if(COMPILER_SUPPORTS_CXX11)
	    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
	elseif(COMPILER_SUPPORTS_CXX0X)
	    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
	else()
		message(FATAL_ERROR "This version of GROMACS requires C++11-compatible standard library. Please use a newer compiler, or a newer standard library, or use the GROMACS 5.1.x release. See the installation guide for details.")
	endif()

endif()
##############################################################

if( HAVE_GROMACS45 )
	add_executable(do_x3dna GMX46/do_x3dna.c ExtractData.h ExtractData.c)
	target_link_libraries(do_x3dna ${GROMACS_LIBRARIES} )
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
	install(TARGETS do_x3dna DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin)
endif()

if( HAVE_GROMACS50 )
	add_executable(do_x3dna GMX50/do_x3dna.c ExtractData.h ExtractData.c)
	target_link_libraries(do_x3dna ${GROMACS_LIBRARIES} )
	SET_TARGET_PROPERTIES(do_x3dna PROPERTIES LINKER_LANGUAGE CXX)
	install(TARGETS do_x3dna DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin)
endif()

if( HAVE_GROMACS51 )
	add_executable(do_x3dna GMX51/do_x3dna.c ExtractData.h ExtractData.c)
	target_link_libraries(do_x3dna ${GROMACS_LIBRARIES} )
	SET_TARGET_PROPERTIES(do_x3dna PROPERTIES LINKER_LANGUAGE CXX)
	install(TARGETS do_x3dna DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin)
endif()

if( HAVE_GROMACS2016 )
	add_executable(do_x3dna GMX2016/do_x3dna.cpp ExtractData.h ExtractData.c)
	target_link_libraries(do_x3dna ${GROMACS_LIBRARIES} )
	SET_TARGET_PROPERTIES(do_x3dna PROPERTIES LINKER_LANGUAGE CXX)
	install(TARGETS do_x3dna DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin)
endif()
