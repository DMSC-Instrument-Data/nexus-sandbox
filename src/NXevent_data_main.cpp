#include <stdio.h>
#include <xmmintrin.h>
#include <mpi.h>
#include <algorithm>
#include <string>
#include <future>
#include <thread>
#include <tuple>

#include "NXevent_data_loader.h"
#include "timer.h"

constexpr size_t bank_size = 1000;

struct Event {
  int32_t index; // global spectrum index
  int32_t tof;
  int64_t pulse_time;
};

struct LoadRange {
  int bank;
  hsize_t start;
  hsize_t count;
};

// In-place translation from event_id to GlobalSpectrumIndex. Note that
// GlobalSpectrumIndex is actually 64 bit, but if the event_id is 32 bit it is
// guaranteed that the index also fits into 32 bit (unless we have scanning,
// which we do not support here for now). Thus we keep things in 32 bit for
// performance.
void event_id_to_GlobalSpectrumIndex(const int32_t event_id_offset,
                                      std::vector<int32_t> &data) {
  // For now we assume that the event_id is contiguous within a bank, so
  // translation is just a fixed offset.
  std::for_each(data.begin(), data.end(),
                [&](int32_t &i) { i -= event_id_offset; });
}

void build_event_vector(std::vector<Event> &events, const LoadRange &range,
                        const std::vector<int32_t> &event_index,
                        const std::vector<int64_t> &event_time_zero,
                        const std::vector<int32_t> &event_global_spectrum_index,
                        const std::vector<int32_t> &event_time_offset) {
  size_t start_pulse = 0;
  size_t end_pulse = 0;
  for (size_t pulse = 0; pulse < event_index.size(); ++pulse) {
    size_t count =
        (pulse != event_index.size() - 1 ? event_index[pulse + 1]
                                         : event_time_offset.size()) -
        event_index[pulse];
    if (range.start >= event_index[pulse] &&
        range.start < event_index[pulse] + count)
      start_pulse = pulse;
    if (range.start + range.count > event_index[pulse] &&
        range.start + range.count <= event_index[pulse] + count)
      end_pulse = pulse + 1;
  }
  // event_index should have offset range.start
  // count 
  events.resize(event_time_offset.size());
  for (size_t pulse = start_pulse; pulse < end_pulse; ++pulse) {
    size_t start =
        std::max(range.start, static_cast<hsize_t>(event_index[pulse])) -
        range.start;
    size_t end =
        std::min(range.start + range.count,
                 static_cast<hsize_t>(pulse != event_index.size() - 1
                                          ? event_index[pulse + 1]
                                          : event_time_offset.size())) -
        range.start;

    for (size_t i = start; i < end; ++i) {
      events[i] = {event_global_spectrum_index[i], event_time_offset[i],
                   event_time_zero[pulse]};
    }
  }
}

template <class T>
void split_for_ranks(std::vector<std::vector<T>> &rank_data, const int nrank,
                     const std::vector<T> &data) {
  for (auto &item : rank_data)
    item.clear();
  for (const auto &item : data) {
    int rank = item.index % nrank;
    rank_data[rank].push_back(item);
  }
}

template <class T>
void redistribute_data(std::vector<T> &result,
                       const std::vector<std::vector<T>> &data) {
  std::vector<size_t> sizes(data.size());
  std::transform(data.cbegin(), data.cend(), sizes.begin(),
                 [](const std::vector<T> &vec) { return vec.size(); });
  std::vector<size_t> rec_sizes(data.size());
  MPI_Alltoall(sizes.data(), 1, MPI_UNSIGNED_LONG, rec_sizes.data(), 1,
               MPI_UNSIGNED_LONG, MPI_COMM_WORLD);

  size_t total_size = std::accumulate(rec_sizes.begin(), rec_sizes.end(), 0);
  result.resize(total_size);
  size_t offset = 0;
  int this_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &this_rank);
  std::vector<MPI_Request> recv_requests(data.size());
  for (int rank = 0; rank < data.size(); ++rank) {
    int tag = 0;
    // Receive from rank with offset to balance network load.
    int rank2 = (rank + this_rank)%data.size();
    MPI_Irecv(result.data() + offset, rec_sizes[rank2] * sizeof(T), MPI_CHAR,
             rank2, tag, MPI_COMM_WORLD, &recv_requests[rank2]);
    offset += rec_sizes[rank2];
    // TODO
    // 2. do work between sending and receiving (next range, or insert into
    // workspace after first data was received?).
  }

  std::vector<MPI_Request> send_requests(data.size());
  for (int rank = 0; rank < data.size(); ++rank) {
    int rank2 = (rank + this_rank)%data.size();
    const auto &vec = data[rank2];
    int tag = 0;
    MPI_Isend(vec.data(), vec.size() * sizeof(T), MPI_CHAR, rank2, tag,
              MPI_COMM_WORLD, &send_requests[rank2]);
  }

  MPI_Waitall(data.size(), send_requests.data(), MPI_STATUSES_IGNORE);
  MPI_Waitall(data.size(), recv_requests.data(), MPI_STATUSES_IGNORE);
}

struct MantidEvent {
  MantidEvent(double tof, int64_t pulse_time)
      : tof(tof), pulse_time(pulse_time) {}
  double tof;
  int64_t pulse_time;
};

template <class T>
void append_to_workspace(const int nrank,
                         std::vector<std::vector<MantidEvent>> &workspace,
                         const std::vector<T> &events) {
  for (const auto &event : events) {
    auto index = event.index / nrank; // global index to local index
    workspace[index].emplace_back(
        static_cast<double>(event.tof), event.pulse_time);
    _mm_prefetch(&workspace[index].back() + 1, _MM_HINT_T1);
  }
}

size_t readEventCount(H5::H5File &file) {
  size_t n_events{0};
  for (int i = 0; i < 7; ++i) {
    H5::DataSet dataset = file.openDataSet("entry/instrument/events-" +
                                           std::to_string(i) + "/event_id");
    H5::DataSpace dataSpace = dataset.getSpace();
    n_events += dataSpace.getSelectNpoints();
  }
  return static_cast<size_t>(n_events);
}

std::vector<LoadRange> determineLoadRanges(int nrank, int rank, H5::H5File &file) {
  // Work balancing should take into account total number of events and then
  // chunk things up for MPI
  size_t n_events{0};
  std::vector<size_t> bank_sizes;
  for (int i = 0; i < 7; ++i) {
    H5::DataSet dataset = file.openDataSet("entry/instrument/events-" +
                                           std::to_string(i) + "/event_id");
    H5::DataSpace dataSpace = dataset.getSpace();
    bank_sizes.push_back(dataSpace.getSelectNpoints());
    n_events += bank_sizes.back();
  }

  std::vector<size_t> bank_ends(bank_sizes);
  for (size_t i = 1; i < bank_ends.size(); ++i)
    bank_ends[i] += bank_ends[i - 1];
  std::vector<size_t> bank_begins(bank_ends);
  for (size_t i = 0; i < bank_begins.size(); ++i)
    bank_begins[i] -= bank_sizes[i];

  size_t offset = n_events / nrank;
  size_t start = offset * rank;
  size_t last = offset * (rank + 1) - 1;
  if (rank == nrank - 1)
    last += n_events - (offset * nrank);

  auto it = std::lower_bound(bank_ends.begin(),
                                   bank_ends.end(), start);
  auto last_it = std::upper_bound(bank_ends.begin(), bank_ends.end(), last);
  std::vector<LoadRange> ranges;
  while(it <= last_it) {
    int bank =
        static_cast<int>(std::distance(bank_ends.begin(), it));
    hsize_t range_start =
        std::max(start, bank_begins[bank]) - bank_begins[bank];

    hsize_t range_end = std::min(last + 1, bank_ends[bank]) - bank_begins[bank];
    ranges.push_back(LoadRange{bank, range_start, range_end - range_start});
    ++it;
  }
  return ranges;
}

std::vector<LoadRange> determineChunkedLoadRanges(int nrank, int rank, H5::H5File &file) {
  size_t chunk_size = 1024*1024; // element count in chunk
  // Work balancing should take into account total number of events and then
  // chunk things up for MPI
  size_t n_events{0};
  std::vector<size_t> bank_sizes;
  for (int i = 0; i < 7; ++i) {
    H5::DataSet dataset = file.openDataSet("entry/instrument/events-" +
                                           std::to_string(i) + "/event_id");
    H5::DataSpace dataSpace = dataset.getSpace();
    bank_sizes.push_back(dataSpace.getSelectNpoints());
    n_events += bank_sizes.back();
  }

  size_t chunk = 0;
  std::vector<LoadRange> ranges;
  for (int bank = 0; bank < 7; ++bank) {
    size_t current = 0;
    while (current < bank_sizes[bank]) {
      if (chunk % nrank == rank) {
        hsize_t count =
            std::min(current + chunk_size, bank_sizes[bank]) - current;
        ranges.push_back(LoadRange{bank, current, count});
      }
      current += chunk_size;
      chunk++;
    }
  }
  while (chunk % nrank != 0) {
    if (chunk % nrank == rank)
        ranges.push_back(LoadRange{0, 0, 0});
    chunk++;
  }
  return ranges;
}

std::tuple<std::vector<int32_t>, std::vector<int64_t>, std::vector<int32_t>,
           std::vector<int32_t>>
load(H5::H5File &file, const LoadRange &range) {
  NXEventDataLoader loader(file, "entry/instrument/events-" +
                                     std::to_string(range.bank));
  // TODO avoid reallocating many times. Reuse buffer (double buffering when
  // threaded?)
  return std::make_tuple(loader.eventIndex(), loader.eventTimeZero(),
          loader.readEventID(range.start, range.count),
          loader.readEventTimeOffset(range.start, range.count));
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int rank;
  int nrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nrank);

  H5::H5File file(argv[1], H5F_ACC_RDONLY);
  size_t n_events = 0;

  //for(const auto &range : ranges)
  //  printf("%d %d %llu %llu\n", rank, range.bank, range.start, range.count);

  int num_loaded_banks = 7;
  std::vector<std::vector<MantidEvent>> workspace(
      (num_loaded_banks * bank_size) / nrank + 1);
  std::vector<Event> events;
  std::vector<std::vector<Event>> events_for_ranks(nrank);
  std::vector<Event> &events_for_this_rank = events;
  std::vector<double> deltas(6, 0.0);

  std::future<std::tuple<std::vector<int32_t>, std::vector<int64_t>,
                         std::vector<int32_t>, std::vector<int32_t>>> future[2];
  const auto ranges = determineChunkedLoadRanges(nrank, rank, file);

  int load_index = 0;
  future[load_index] =
      std::async(std::launch::async, load, std::ref(file), ranges[0]);
  load_index = (load_index + 1) % 2;
  for (size_t range_index = 0; range_index < ranges.size(); ++range_index) {
    const auto range = ranges[range_index];
    Timer timer;

    if (range_index != ranges.size() - 1)
      future[load_index] = std::async(std::launch::async, load, std::ref(file),
                                      ranges[range_index + 1]);
    load_index = (load_index + 1) % 2;

    future[load_index].wait();
    auto result = future[load_index].get();
    auto event_index = std::get<0>(result);
    auto event_time_zero = std::get<1>(result);
    auto event_id = std::get<2>(result);
    auto event_time_offset = std::get<3>(result);
    timer.checkpoint();

    event_id_to_GlobalSpectrumIndex(10 * bank_size * range.bank, event_id);
    timer.checkpoint();
    // event_id is now actually event_global_spectrum_index
    build_event_vector(events, range, event_index, event_time_zero,
                                            event_id, event_time_offset);
    timer.checkpoint();
    split_for_ranks(events_for_ranks , nrank, events);
    timer.checkpoint();
    redistribute_data(events_for_this_rank, events_for_ranks);
    timer.checkpoint();
    append_to_workspace(nrank, workspace, events_for_this_rank);
    timer.checkpoint();
    const auto seconds = timer.deltas();
    for(size_t i=0; i<deltas.size(); ++i)
      deltas[i] += seconds[i];
  }

  size_t n_event_in_local_workspace{0};
  size_t n_event_in_workspace{0};
  for (const auto &item : workspace)
    n_event_in_local_workspace += item.size();
  MPI_Reduce(&n_event_in_local_workspace, &n_event_in_workspace, 1,
             MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    printf("HDF5-load id-to-index event-vec split-ranks MPI make-worksapce\n");
    double sum{0.0};
    for (const auto seconds : deltas) {
      sum += seconds;
      printf("%.3lf ", seconds);
    }
    auto n_event = readEventCount(file);
    printf("total %.2lf | %.1e events/s\n", sum, n_event/sum);
    if (n_event != n_event_in_workspace)
      throw std::runtime_error("lost events\n");
  }

  MPI_Finalize();
  return 0;
}
