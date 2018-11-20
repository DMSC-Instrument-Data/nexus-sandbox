import argparse
import h5py
import numpy as np
from shutil import copyfile
from datetime import datetime

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i", "--input-filename", type=str, help='Input file to convert.')
parser.add_argument("-r", "--reference-filename", type=str, help='Input file used as a reference for the file structure.')
parser.add_argument("-o", "--output-filename", type=str, help='Output filename.')
args = parser.parse_args()

input_config = dict(root='entry-01', event_data=['Delayline_events'], Nx=16, Ny=32)
# Used PG3_4866_event.nxs for testing
output_config = dict(root='entry', event_data=['bank102_events'], event_id_offset=[128000])

def position_to_index(pos, count):
    uint_max = 2**16-1
    # What if count (Nx or Ny) does not divide uint_max?
    return np.floor_divide(pos, (uint_max//count))

def convert_id(event_id, id_offset):
    x = (event_id[:] >> 16) & 0xffff
    y = event_id[:] & 0xffff
    Nx = input_config['Nx']
    Ny = input_config['Ny']
    # Mantid requires 32 bit unsigned, so this should be correct dtype already.
    # Need offset here unless the banks event ids start at zero (Mantid
    # will simply discard events that do not correspond to IDF).
    return id_offset + position_to_index(x, Nx) + Nx * position_to_index(y, Ny)

def make_index(absolute_times, time_zero):
    # Mantid uses 64 bit unsigned int for event_index
    index = np.zeros(shape=time_zero.shape, dtype=np.uint64)
    cur = 1
    for i, t in enumerate(absolute_times):
        # Assuming that absolute_times and time_zero has the same offset.
        if t >= time_zero[cur]:
            index[cur] = i
            cur += 1
            if cur == len(index):
                break
    return index

def to_seconds(nanoseconds):
    return (nanoseconds[:]/1e9).astype(dtype=np.float64)

def convert_time(absolute_times):
    absolute_times = to_seconds(absolute_times)
    # Offset stored as attribute of event_time_zero. What is the base in our files?
    time_zero_offset = datetime(year=2018, month=1, day=1).isoformat()
    # TODO Currently I do not know where to get this from (parse TDC?), creating dummy data for testing.
    # Mantid uses 64 bit event_time_zero with unit second.
    time_zero = np.arange(0.0, 100.0, 0.071, dtype=np.float64)
    # Mantid uses 32 bit event_time_offset with unit microsecond.
    time_offset = np.random.uniform(0.0, 71000.0, size=absolute_times.shape).astype(dtype=np.float32)
    index = make_index(absolute_times, time_zero)
    return time_zero_offset, time_zero, time_offset, index

if __name__ == '__main__':
    copyfile(args.reference_filename, args.output_filename)

    with h5py.File(args.input_filename, 'r') as input_file, h5py.File(args.output_filename, 'a') as output_file:
        # Delete unused banks from reference file
        # Note: This does not actually make the file smaller.
        for bank, event_data in output_file[output_config['root']].items():
            if bank.endswith('_events') and not bank in output_config['event_data']:
                del output_file[output_config['root']][bank]

        for src, dst, id_offset in zip(input_config['event_data'], output_config['event_data'], output_config['event_id_offset']):
            src_data = input_file[input_config['root']][src]
            dst_data = output_file[output_config['root']][dst]
            print('Loading event data from group {}, writing into {}'.format(src, dst))

            # event_index is just a monotonically increasing integer, not needed.
            # event_index = data['event_index']

            # event_id is 16 bit x pos and 16 bit y pos -> convert
            event_id = convert_id(src_data['event_id'], id_offset)

            # event_time_offset is zero or amplitude of pulse, not needed.
            # event_time_offset = data['event_time_offset']

            # event_time_zero is the absolute time stamp of each event.
            absolute_times = src_data['event_time_zero']
            time_zero_offset, time_zero, time_offset, index = convert_time(absolute_times)

            dst_data['event_id'].resize(len(event_id), axis = 0)
            dst_data['event_id'][:] = event_id

            dst_data['event_time_zero'].resize(len(time_zero), axis = 0)
            dst_data['event_time_zero'][:] = time_zero
            dst_data['event_time_zero'].attrs['offset'] = time_zero_offset

            dst_data['event_time_offset'].resize(len(time_offset), axis = 0)
            dst_data['event_time_offset'][:] = time_offset

            dst_data['event_index'].resize(len(index), axis = 0)
            dst_data['event_index'][:] = index
