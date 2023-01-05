cmake_minimum_required(VERSION 3.2)
#policy CMP0076 - target_sources source files are relative to file where target_sources is run
cmake_policy (SET CMP0076 NEW)

# MiMA claims to only compile with ifort / icc currently.
set ( CMAKE_Fortran_COMPILER "ifort" )
set ( CMAKE_C_COMPILER "icc" )
project(MiMA Fortran C)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "" FORCE)
endif()

# Generate the main mima.x executable with dirs, libs, and opts
add_executable ( mima.x )
target_include_directories( mima.x PUBLIC src/shared/include src/shared/mpp/include )
add_library( mima_c ) # The C parts of MiMA, so we can apply different options for them.
target_compile_definitions( mima_c PRIVATE __IFC )
target_compile_definitions( mima.x PRIVATE use_libMPI use_netCDF gFortran ) # gFortran appears to be unused

# Also generate the postprocessing executable
add_executable ( mppnccombine postprocessing/mppnccombine.c )

#Add cmake directory to the environment module variable
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Set to install in bin directory as per current MiMA behaviour
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/bin" CACHE PATH "..." FORCE)
endif()

# Find MPI, OpenMP, and python libraries and link
find_package (MPI REQUIRED)
target_link_libraries( mima.x PRIVATE MPI::MPI_Fortran )

find_package (Python REQUIRED COMPONENTS Development)
target_link_libraries( mima.x PRIVATE Python::Python )

find_package (OpenMP REQUIRED COMPONENTS Fortran)
target_link_libraries( mima.x PRIVATE OpenMP::OpenMP_Fortran )

# Find the NetCDF installations and set the relevant variables for compilation
# Then link to executables
# Requires more legwork as NetCDF not provided by default
find_package(PkgConfig)
pkg_search_module(NETCDF_FORTRAN netcdf-fortran)
if (NETCDF_FORTRAN_FOUND)
    set(NETCDF_LIBRARIES "${NETCDF_FORTRAN_LDFLAGS}")
    set(NETCDF_INCLUDES "${NETCDF_FORTRAN_INCLUDE_DIRS}")
else()
    set(NETCDF_F90 "YES")
    find_package(NetCDF REQUIRED)
endif()
pkg_search_module(NETCDF_C netcdf)
if (NETCDF_C_FOUND)
    set(NETCDF_LIBRARIES "${NETCDF_LIBRARIES} ${NETCDF_C_LDFLAGS}")
    set(NETCDF_INCLUDES "${NETCDF_LIBRARIES} ${NETCDF_C_INCLUDE_DIRS}")
endif()

target_link_libraries( mima.x PRIVATE mima_c ${NETCDF_LIBRARIES} )
target_include_directories( mima.x PRIVATE ${NETCDF_INCLUDES} )
target_link_libraries( mppnccombine PRIVATE ${NETCDF_LIBRARIES} )
target_include_directories( mppnccombine PRIVATE ${NETCDF_INCLUDES} )

# Add various subdirectories with long lists of source files
add_subdirectory( src/coupler )
add_subdirectory( src/atmos_coupled )
add_subdirectory( src/atmos_param )

set_source_files_properties (
# The following files do nothing but assign very large arrays.
# For some reason when compiling with ifort and optimisation
# the compilation will take a very long time (10s of minutes).
# Since the code doesn't actually *do* anything there's no
# need to waste time having the compiler apply probably
# meaningless optimisation.  So we disable optimisation for these
# files.
src/atmos_param/rrtm_radiation/rrtmg_sw/gcm_model/src/rrtmg_sw_k_g.f90
src/atmos_param/rrtm_radiation/rrtmg_lw/gcm_model/src/rrtmg_lw_k_g.f90
PROPERTIES
COMPILE_FLAGS -O0
)
add_subdirectory( src/atmos_shared )
add_subdirectory( src/atmos_spectral )
add_subdirectory( src/ice_param )
# include/fms_platform.h
add_subdirectory( src/shared )

# Set coompile options for executable
target_compile_options( mima.x PRIVATE
-fpp
-safe-cray-ptr
-ftz
-assume byterecl
-i4
-r8
)
