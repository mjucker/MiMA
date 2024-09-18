cmake_minimum_required(VERSION 3.16)
#policy CMP0076 - target_sources source files are relative to file where target_sources is run
cmake_policy (SET CMP0076 NEW)
project(MiMA Fortran C)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "" FORCE)
endif()

# Generate the main mima.x executable with dirs, libs, and opts
add_executable ( mima.x )
target_include_directories( mima.x PUBLIC src/shared/include src/shared/mpp/include )
add_library( mima_c ) # The C parts of MiMA, so we can apply different options for them.
target_compile_definitions( mima.x PRIVATE use_libMPI use_netCDF gFortran ) # gFortran appears to be unused

# Also generate the postprocessing executable
add_executable ( mppnccombine postprocessing/mppnccombine.c )

#Add cmake directory to the environment module variable
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Set to install in bin directory as per current MiMA behaviour
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/exec" CACHE PATH "..." FORCE)
endif()

# Find MPI, OpenMP, and python libraries and link
find_package (MPI REQUIRED)
target_link_libraries( mima.x PRIVATE MPI::MPI_Fortran )

find_package (Python REQUIRED COMPONENTS Development)
target_link_libraries( mima.x PRIVATE Python::Python )

find_package (OpenMP REQUIRED COMPONENTS Fortran)
target_link_libraries( mima.x PRIVATE OpenMP::OpenMP_Fortran )

find_package(NetCDF REQUIRED COMPONENTS C Fortran)
target_link_libraries(mima.x PRIVATE NetCDF::NetCDF_Fortran NetCDF::NetCDF_C)
target_link_libraries(mppnccombine PRIVATE NetCDF::NetCDF_C)

if(NOT CMAKE_Fortran_COMPILER)
    set(CMAKE_Fortran_COMPILER "gfortran")
endif()
if(MPI_Fortran_COMPILER)
    set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
    message(STATUS "Using MPI Fortran compiler: ${MPI_Fortran_COMPILER}")   
endif()


if(APPLE)
    message(STATUS "Configuring for MacOS")
    if(NOT CMAKE_C_COMPILER)
        set(CMAKE_C_COMPILER "clang")
    endif()
    set(CMAKE_OSX_ARCHITECTURES "arm64" CACHE STRING "Build architectures for MacOS")
endif()

if(NOT CMAKE_C_COMPILER)
    set(CMAKE_C_COMPILER "gcc")
endif()


if(CMAKE_Fortran_COMPILER_ID  MATCHES "Intel|IntelLLVM")
    message(STATUS "Configuring for Intel Fortran")
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -i4 -r8 -fpp -fPIC")
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O3 -xHost")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    message(STATUS "Configuring for GNU Fortran")
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -fallow-argument-mismatch -fbacktrace -fdefault-real-8 -fdefault-double-8 -fcray-pointer -fallow-invalid-boz -fno-range-check -ffree-line-length-none")
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O3 -march=native")
endif()

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
add_subdirectory( src/shared )

target_link_libraries(mima.x PRIVATE mima_c)
install(TARGETS mima.x mppnccombine RUNTIME DESTINATION .) 
install(FILES input/diag_table DESTINATION . )
install(FILES input/field_table DESTINATION .)
install(DIRECTORY DESTINATION RESTART)
install(FILES input/input.nml DESTINATION . )  
install(DIRECTORY input/INPUT/ DESTINATION .) 