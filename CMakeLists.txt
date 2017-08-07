cmake_minimum_required ( VERSION 3.5 )
project ( nexus-sandbox )

find_package ( HDF5 COMPONENTS CXX HL REQUIRED )

add_subdirectory ( src )
