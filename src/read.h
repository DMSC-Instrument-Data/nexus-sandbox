#ifndef READ_H
#define READ_H

#include <vector>
#include <H5Cpp.h>

template <class T>
std::vector<T> read(const H5::H5File &file, const std::string &dataSetName) {
  H5::DataSet dataset = file.openDataSet(dataSetName);
  H5::DataType dataType = dataset.getDataType();
  H5::DataSpace dataSpace = dataset.getSpace();
  std::vector<T> result;
  result.resize(dataSpace.getSelectNpoints());
  dataset.read(result.data(), dataType);
  return result;
}

template <class T>
void read(std::vector<T> &result, const H5::H5File &file,
          const std::string &dataSetName, hsize_t start, hsize_t count) {
  H5::DataSet dataset = file.openDataSet(dataSetName);
  H5::DataType dataType = dataset.getDataType();
  H5::DataSpace dataSpace = dataset.getSpace();
  if (static_cast<int64_t>(dataSpace.getSelectNpoints()) -
          static_cast<int64_t>(start) <=
      0)
    return;
  count = std::min(count, dataSpace.getSelectNpoints() - start);
  dataSpace.selectHyperslab(H5S_SELECT_SET, &count, &start);
  result.resize(dataSpace.getSelectNpoints());
  H5::DataSpace memSpace(1, &count);
  dataset.read(result.data(), dataType, memSpace, dataSpace);
}

template <class T>
std::vector<T> read(const H5::H5File &file, const std::string &dataSetName,
                    hsize_t start, hsize_t count) {
  std::vector<T> result;
  read(result, file, dataSetName, start, count);
  return result;
}

#endif // READ_H
