#include <stdio.h>
#include <chrono>
#include <vector>

#include <boost/mpi.hpp>
#include <H5Cpp.h>

#include "read.h"

int main(int argc, char **argv) {
  boost::mpi::environment env;
  boost::mpi::communicator comm;

  // /home/simon/mantid/nexus/load-performance/sample-files/no-compression.nxs
  // /home/simon/mantid/nexus/load-performance/sample-files/gzip-level6.nxs
  H5::H5File file(argv[1], H5F_ACC_RDONLY);

  //std::vector<int> banks = {22,  23,  24,  42,  43,  44,  62,  63,
  //                          64,  82,  83,  84,  102, 103, 104, 105,
  //                          106, 123, 124, 143, 144, 164, 184};

  std::vector<int> banks = {28};

  std::string prefix("entry/bank");
  std::string suffix("_events/");
  std::string name = prefix + std::to_string(banks[0]) + suffix;
  //std::string name = prefix + std::to_string(banks[comm.rank()]) + suffix;

  const auto start = std::chrono::system_clock::now();
  size_t size = 0;

  comm.barrier();

  {
  hsize_t count = 2048*1024;
  hsize_t start = count * comm.rank();
  while (true) {
    const auto data = read<int32_t>(file, name + "event_id", start, count);
    if(data.empty())
      break;
    start += count * comm.size();
    size += data.size() * sizeof(int32_t);
  }
  }
  /*
  auto id = read<int32_t>(file, name + "event_id");
  auto tof = read<float>(file, name + "event_time_offset");
  {
  hsize_t start = 0;
  hsize_t count = 2048*1;
  while (read<int32_t>(file, name + "event_id", start, count).size() == count) {
    start += count;
  }
  }
  */
  comm.barrier();

  const auto end = std::chrono::system_clock::now();

  std::chrono::duration<double> elapsed = end - start;
  double seconds = (double)std::chrono::duration_cast<std::chrono::nanoseconds>(
                       elapsed).count() /
                   (1000000000);

  size_t sum;
  boost::mpi::reduce(comm, size, sum, std::plus<size_t>(), 0);
  if (comm.rank() == 0)
    printf("%f bandwidth %f MB/s size %lu\n", seconds,
           static_cast<double>(sum) / seconds / (1024 * 1024), sum);
}
