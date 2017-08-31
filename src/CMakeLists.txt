add_executable ( load_parallel load_parallel.cpp )
target_include_directories ( load_parallel PRIVATE ${HDF5_INCLUDE_DIRS} ${MPI_INCLUDE_PATH} )
target_link_libraries ( load_parallel ${HDF5_CXX_LIBRARIES} ${MPI_LIBRARIES} Boost::mpi )

add_executable ( load_NXevent_data NXevent_data_main.cpp NXevent_data_loader.cpp )
target_include_directories ( load_NXevent_data PRIVATE ${HDF5_INCLUDE_DIRS} ${MPI_INCLUDE_PATH} )
target_link_libraries ( load_NXevent_data ${HDF5_CXX_LIBRARIES} ${MPI_LIBRARIES} Boost::mpi )
