#include <stdio.h>
#include <vector>

#include <mpi.h>
#include <H5Cpp.h>

#include "read.h"
#include "timer.h"

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int rank;
  int world_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // /home/simon/mantid/nexus/load-performance/sample-files/no-compression.nxs
  // /home/simon/mantid/nexus/load-performance/sample-files/gzip-level6.nxs
  // /mnt/extra/simon/neutron-data/VIS_20118.nxs.h5
  H5::H5File file(argv[1], H5F_ACC_RDONLY);

  //std::vector<int> banks = {22,  23,  24,  42,  43,  44,  62,  63,
  //                          64,  82,  83,  84,  102, 103, 104, 105,
  //                          106, 123, 124, 143, 144, 164, 184};

  std::vector<int> banks = {1};

  std::string prefix("entry/bank");
  std::string suffix("_events/");
  std::string name = prefix + std::to_string(banks[0]) + suffix;
  //std::string name = prefix + std::to_string(banks[comm.rank()]) + suffix;

  Timer timer;
  size_t size = 0;

  MPI_Barrier(MPI_COMM_WORLD);

  {
  hsize_t count = 2048*1024;
  hsize_t start = count * rank;
  while (true) {
    const auto data = read<int32_t>(file, name + "event_id", start, count);
    if(data.empty())
      break;
    start += count * world_size;
    size += data.size() * sizeof(int32_t);
  }
  /*
    const auto data = read<int32_t>(file, name + "event_id");
    size += data.size() * sizeof(int32_t);
  */
  }

  MPI_Barrier(MPI_COMM_WORLD);

  double seconds = timer.checkpoint();

  size_t sum;
  MPI_Reduce(&size, &sum, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
  if (rank == 0)
    printf("%f bandwidth %f MB/s size %lu\n", seconds,
           static_cast<double>(sum) / seconds / (1024 * 1024), sum);

  MPI_Finalize();
}
