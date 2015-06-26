# ============== Eclipse CDT Error Handling======
# IF(CMAKE_COMPILER_IS_GNUCC)
#  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fmessage-length=0")
# ENDIF(CMAKE_COMPILER_IS_GNUCC)
# IF(CMAKE_COMPILER_IS_GNUCXX)
#  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmessage-length=0")
# ENDIF(CMAKE_COMPILER_IS_GNUCXX)
# set(CMAKE_VERBOSE_MAKEFILE ON)
# ===============================

# Checking GROMACS PAth
if (NOT DEFINED GMX_PATH AND (NOT DEFINED GMX_LIB OR NOT DEFINED GMX_INCLUDE) )
message(FATAL_ERROR "\nGromacs path is not defined.....\nUSE:\n -DGMX_PATH=/path/to/gromacs \nOR\n -DGMX_LIB=/path/to/gromacs/lib -DGMX_INCLUDE=/path/to/gromacs/include")
endif()


if(DEFINED GMX_PATH)
	set(GMX_INCLUDE "${GMX_PATH}/include")
	set(GMX_LIB "${GMX_PATH}/lib")
endif()

include_directories("${GMX_INCLUDE}/gromacs")

# Searching FFTW library
find_library(LIBFFTW fftw3f HINTS ${FFTW_LIB})
if(DEFINED LIBFFTW-NOTFOUND)
message(FATAL_ERROR "\nFFTW library file libfftw3f.so or libfftw3f.a not found at standard location\nUSE:\n -DFFTW_LIB=/path/to/fftw3/lib")
endif()

# Searching GROMACS and other Library files
find_library(LIBM m)
find_library(LIBGMX gmx HINTS ${GMX_LIB})
find_library(LIBGMXANA gmxana HINTS ${GMX_LIB})
find_library(LIBMD md HINTS ${GMX_LIB})
find_library(LIBPTHREAD pthread)
find_library(LIBDL dl)

add_executable(do_x3dna do_x3dna.c ExtractData.h ExtractData.c)
target_link_libraries(do_x3dna ${LIBGMX} ${LIBGMXANA} ${LIBMD} ${LIBFFTW} ${LIBM} ${LIBDL} ${LIBPTHREAD} )
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
install(TARGETS do_x3dna DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin)
