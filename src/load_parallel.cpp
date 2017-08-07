#include <stdio.h>
#include <chrono>
#include <vector>
#include <H5Cpp.h>

template <class T>
std::vector<T> read(const H5::H5File &file, const std::string &dataSetName) {
  H5::DataSet dataset = file.openDataSet(dataSetName);
  H5::DataType dataType = dataset.getDataType();
  H5::DataSpace dataSpace = dataset.getSpace();
  std::vector<float> result;
  result.resize(dataSpace.getSelectNpoints());
  dataset.read(result.data(), dataType, dataSpace);
  return result;
}

int main() {
  H5::H5File file(
      "/home/simon/mantid/nexus/load-performance/sample-files/gzip-level6.nxs",
      H5F_ACC_RDONLY);

  const auto start = std::chrono::system_clock::now();

  auto result1 = read<float>(file, "entry/bank22_events/event_time_offset");
  auto result2 = read<float>(file, "entry/bank23_events/event_time_offset");

  const auto end = std::chrono::system_clock::now();

  std::chrono::duration<double> elapsed = end - start;
  double seconds = (double)std::chrono::duration_cast<std::chrono::nanoseconds>(
                       elapsed).count() /
                   (1000000000);

  printf("%f bandwidth %f MB/s\n", seconds,
         static_cast<double>((result1.size() + result2.size()) * sizeof(float)) / seconds /
             (1024 * 1024));
}
