#include <stdio.h>
#include <string>
#include "NXevent_data_loader.h"
#include "timer.h"

int main(int argc, char **argv) {
  H5::H5File file(argv[1], H5F_ACC_RDONLY);
  size_t n_events = 0;
  Timer timer;
  // Work balancing should take into account total number of events and then
  // chunk things up for MPI
  for (int i = 0; i < 7; ++i) {
    H5::DataSet dataset = file.openDataSet("entry/instrument/events-" +
                                           std::to_string(i) + "/event_id");
    H5::DataSpace dataSpace = dataset.getSpace();
    n_events += dataSpace.getSelectNpoints();
  }
  printf("counted %lu events total %lf seconds\n", n_events, timer.checkpoint());
  NXEventDataLoader loader(file, "entry/instrument/events-0");
  // TODO avoid reallocating many times. Reuse buffer (double buffering when
  // threaded?)
  auto ids = loader.readEventID(0, 10);
  return 0;
}
