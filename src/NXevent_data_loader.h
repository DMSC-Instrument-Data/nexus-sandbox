#ifndef NXEVENT_DATA_LOADER_H
#define NXEVENT_DATA_LOADER_H

#include <string>
#include <vector>

#include <H5Cpp.h>

class NXEventDataLoader {
public:
  NXEventDataLoader(const H5::H5File &file, const std::string &nxEventDataPath);

  size_t numberOfEvents() const;
  const std::vector<int32_t> &eventIndex() const;
  const std::vector<int64_t> &eventTimeZero() const;
  std::vector<int32_t> readEventID(hsize_t start, hsize_t count) const;
  std::vector<int32_t> readEventTimeOffset(hsize_t start, hsize_t count) const;

private:
  const H5::H5File &m_file;
  const std::string m_NXEventDataPath;
  const std::string m_event_id_path;
  const std::string m_event_time_offset_path;
  std::vector<int32_t> m_event_index;
  std::vector<int64_t> m_event_time_zero;
};

#endif // NXEVENT_DATA_LOADER_H
