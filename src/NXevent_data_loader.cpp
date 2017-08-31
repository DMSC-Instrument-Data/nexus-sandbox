#include "NXevent_data_loader.h"
#include "read.h"

NXEventDataLoader::NXEventDataLoader(const H5::H5File &file,
                                     const std::string &nxEventDataPath)
    : m_file(file), m_NXEventDataPath(nxEventDataPath),
      m_event_id_path(m_NXEventDataPath + "/event_id"),
      m_event_time_offset_path(m_NXEventDataPath + "/event_time_offset") {
  m_event_index = read<int32_t>(file, m_NXEventDataPath + "/event_index");
  m_event_time_zero = read<int64_t>(file, m_NXEventDataPath + "/event_time_zero");
}

const std::vector<int32_t> &NXEventDataLoader::eventIndex() const {
  return m_event_index;
}
const std::vector<int64_t> &NXEventDataLoader::eventTimeZero() const {
  return m_event_time_zero;
}

std::vector<int32_t> NXEventDataLoader::readEventID(hsize_t start,
                                                    hsize_t count) const {
  return read<int32_t>(m_file, m_event_id_path, start, count);
}

std::vector<int32_t>
NXEventDataLoader::readEventTimeOffset(hsize_t start, hsize_t count) const {
  return read<int32_t>(m_file, m_event_time_offset_path, start, count);
}
