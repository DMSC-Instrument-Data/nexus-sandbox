#include <stdio.h>
#include <mpi.h>
#include <algorithm>
#include <string>

#include "NXevent_data_loader.h"
#include "timer.h"

struct Event {
  int32_t index; // global spectrum index
  int32_t tof;
  int64_t pulse_time;
};

// In-place translation from event_id to GlobalSpectrumIndex. Note that
// GlobalSpectrumIndex is actually 64 bit, but if the event_id is 32 bit it is
// guaranteed that the index also fits into 32 bit (unless we have scanning,
// which we do not support here for now). Thus we keep things in 32 bit for
// performance.
void event_id_to_GlobalSpectrumIndex(const int32_t event_id_offset,
                                      std::vector<int32_t> data) {
  // For now we assume that the event_id is contiguous within a bank, so
  // translation is just a fixed offset.
  std::for_each(data.begin(), data.end(),
                [&](int32_t &i) { i -= event_id_offset; });
}

std::vector<Event>
build_event_vector(const std::vector<int32_t> &event_index,
                   const std::vector<int64_t> &event_time_zero,
                   const std::vector<int32_t> &event_global_spectrum_index,
                   const std::vector<int32_t> &event_time_offset) {
  std::vector<Event> events(event_time_offset.size());
  for (size_t pulse = 0; pulse < event_index.size(); ++pulse) {
    size_t base = event_index[pulse];
    size_t count = (pulse != event_index.size() - 1
                        ? event_index[pulse + 1]
                        : (event_time_offset.size()) - event_index[pulse]);
    for (size_t offset = 0; offset < count; ++offset) {
      size_t i = base + offset;
      events[i] = {event_global_spectrum_index[i], event_time_offset[i],
                   event_time_zero[pulse]};
    }
  }
  return events;
}

template <class T>
std::vector<std::vector<T>> split_for_ranks(const int nrank,
                                              std::vector<T> data) {
  std::vector<std::vector<T>> rank_data(nrank);
  for (const auto &item : data) {
    int rank = item.index % nrank;
    rank_data[rank].push_back(item);
  }
  return rank_data;
}

template <class T>
std::vector<T>
redistribute_data(const std::vector<std::vector<T>> &data) {
  std::vector<size_t> sizes(data.size());
  std::transform(data.cbegin(), data.cend(), sizes.begin(),
                 [](const std::vector<T> &vec) { return vec.size(); });
  std::vector<size_t> rec_sizes(data.size());
  MPI_Alltoall(sizes.data(), 1, MPI_UINT64_T, rec_sizes.data(), 1, MPI_UINT64_T,
               MPI_COMM_WORLD);

  std::vector<MPI_Request> send_requests(data.size());
  for (int rank = 0; rank < data.size(); ++rank) {
    const auto &vec = data[rank];
    int tag = 0;
    MPI_Isend(vec.data(), vec.size() * sizeof(T), MPI_CHAR, rank, tag, MPI_COMM_WORLD,
              &send_requests[rank]);
  }

  size_t total_size = std::accumulate(rec_sizes.begin(), rec_sizes.end(), 0);
  std::vector<T> result(total_size);
  size_t offset = 0;
  for (int rank = 0; rank < data.size(); ++rank) {
    int tag = 0;
    MPI_Recv(result.data() + offset, rec_sizes[rank] * sizeof(T), MPI_CHAR,
             rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    offset += rec_sizes[rank];
  }

  MPI_Waitall(data.size(), send_requests.data(), MPI_STATUSES_IGNORE);
  return result;
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int rank;
  int nrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nrank);

  /*
  std::vector<std::vector<int32_t>> data(nrank);
  for(size_t i=0; i<data.size(); ++i)
    data[i].resize(rank, 1000*rank + i);
  const auto &result = redistribute_data(data);
  for(const auto &i:result)
    printf("%d %d\n", rank, i);
    */







  H5::H5File file(argv[1], H5F_ACC_RDONLY);
  size_t n_events = 0;
  Timer timer;
  /*
  // Work balancing should take into account total number of events and then
  // chunk things up for MPI
  for (int i = 0; i < 7; ++i) {
    H5::DataSet dataset = file.openDataSet("entry/instrument/events-" +
                                           std::to_string(i) + "/event_id");
    H5::DataSpace dataSpace = dataset.getSpace();
    n_events += dataSpace.getSelectNpoints();
  }
  printf("counted %lu events total %lf seconds\n", n_events,
         timer.checkpoint());
         */
  NXEventDataLoader loader(file,
                           "entry/instrument/events-" + std::to_string(rank));
  // TODO avoid reallocating many times. Reuse buffer (double buffering when
  // threaded?)
  const size_t n_event = loader.numberOfEvents();
  const auto &event_index = loader.eventIndex();
  const auto &event_time_zero = loader.eventTimeZero();
  const auto &event_id = loader.readEventID(0, n_event);
  const auto &event_time_offset = loader.readEventTimeOffset(0, n_event);

  event_id_to_GlobalSpectrumIndex(100000, event_id);
  // event_id is now actually event_global_spectrum_index
  const auto &events = build_event_vector(event_index, event_time_zero, event_id, event_time_offset);
  const auto &events_for_ranks = split_for_ranks(nrank, events);
  const auto &event_for_this_rank = redistribute_data(events_for_ranks);

  for(size_t i=0; i<10; ++i) {
    const auto event = event_for_this_rank[i];
    printf("%d %d %d %ld\n", rank, event.index, event.tof, event.pulse_time);
  }

  MPI_Finalize();
  return 0;
}
