cmake_minimum_required ( VERSION 3.5 )
project ( nexus-sandbox )

set ( CMAKE_CXX_STANDARD 11 )
set( CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE )

find_package ( HDF5 COMPONENTS CXX HL REQUIRED )
find_package ( MPI REQUIRED )

add_subdirectory ( src )
