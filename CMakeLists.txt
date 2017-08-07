cmake_minimum_required ( VERSION 3.5 )
project ( nexus-sandbox )

set ( CMAKE_CXX_STANDARD 11 )

find_package ( HDF5 COMPONENTS CXX HL REQUIRED )
find_package ( MPI REQUIRED )
find_package ( Boost REQUIRED mpi )

add_subdirectory ( src )
