import argparse
import h5py
import numpy as np
from shutil import copyfile
from datetime import datetime

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i", "--input-filename", type=str, help='Input file to convert.')
parser.add_argument("-m", "--metadata-filename", type=str, help='Input file used as source for chopper TDC information.')
parser.add_argument("-r", "--reference-filename", type=str, help='Input file used as a reference for the file structure.')
parser.add_argument("-o", "--output-filename", type=str, help='Output filename.')
args = parser.parse_args()

# Example usage:
# python3 convert_nonstandard.py -i V20_example_1.nxs -m V20_example_1.nxs -r LOKI_definition.hdf5 -o test.nxs

# event_data is a list of bank names to convert.
# For adc_test6.nxs (from Jonas)
# input_config = dict(root='entry-01', event_data=['Delayline_events'], Nx=150, Ny=150)
# For V20_example_1.nxs (from Matt)
input_config = dict(root='entry/instrument', event_data=['detector_1/raw_event_data', 'detector_1/raw_event_data'], Nx=150, Ny=150)

# HDF5 group names giving a path to the TDC used as pulse times.
metadata_config = dict(NXdisk_chopper='entry/instrument/chopper_1')

# Used LOKI_definition.hdf5 testing
# NXdetector is a list of output bank names (must exist in reference file), length must match that in input_config.
# event_id_offset must have same length as event_data. This is the offset of the detector ID or spectrum number for the given bank.
#output_config = dict(root='entry', event_data=['bank102_events'], event_id_offset=[128000])
# event_id_offset is bank*1000 for LOKI
output_config = dict(root='raw_data_1', NXinstrument='instrument', NXdetector=['detector_1', 'detector_2'], event_id_offset=[0, 1000])

def position_to_index(pos, count):
    uint_max = 2**16-1
    # What if count (Nx or Ny) does not divide uint_max?
    return np.floor_divide(pos, (uint_max//count))

def convert_id(event_id, id_offset):
    # TODO Is the order correct? Is x in the high bits or the low bits?
    x = (event_id[:] >> 16) & 0xffff
    y = event_id[:] & 0xffff
    Nx = input_config['Nx']
    Ny = input_config['Ny']
    # Mantid requires 32 bit unsigned, so this should be correct dtype already.
    # Need offset here unless the banks event ids start at zero (Mantid
    # will simply discard events that do not correspond to IDF).
    return id_offset + position_to_index(x, Nx) + Nx * position_to_index(y, Ny)

def make_index_and_offset(absolute_times, time_zero):
    # Mantid uses 64 bit unsigned int for event_index
    index = np.zeros(shape=time_zero.shape, dtype=np.uint64)
    time_offset = absolute_times.copy()
    cur = 0
    for i, t in enumerate(absolute_times):
        # Assuming that absolute_times and time_zero has the same offset.
        while cur < len(index)-1 and t >= time_zero[cur+1]:
            index[cur] = i
            cur += 1
        time_offset[i] -= time_zero[cur]
    # Mantid uses 32 bit event_time_offset with unit microsecond.
    time_offset = (time_offset*1e6).astype(dtype=np.float32)
    return index, time_offset

def to_seconds(nanoseconds):
    return (nanoseconds[:]/1e9).astype(dtype=np.float64)

def to_iso8601(nanoseconds):
    dt = datetime.fromtimestamp(nanoseconds // 1000000000)
    s = dt.isoformat() + '.' + str(int(nanoseconds % 1000000000)).zfill(9)
    return s

def unix_epoch_to_epics_epoch_offset():
    # TODO: Note that this ignores leap seconds, which is probably an issue.
    return int((datetime(year=1990, month=1, day=1) - datetime(year=1970, month=1, day=1)).total_seconds())

def convert_time(absolute_times):
    # Mantid uses 64 bit event_time_zero with unit second.
    with h5py.File(args.metadata_filename, 'r') as metadata_file:
        time_zero = metadata_file[metadata_config['NXdisk_chopper']]['top_dead_centre/time']
        # Offset to be stored as attribute of event_time_zero.
        time_zero_offset = time_zero[0]
        end_time = to_iso8601(time_zero[-1])
        time_zero -= time_zero_offset
        time_zero = to_seconds(time_zero)
        # TODO In our current test files the absolute times appear to have a different offset.
        # TODO It *might* be since beginning of the day?
        absolute_times = absolute_times[:] + unix_epoch_to_epics_epoch_offset()
        absolute_times -= time_zero_offset
        absolute_times = to_seconds(absolute_times)
        time_zero_offset = to_iso8601(time_zero_offset)
        index, time_offset = make_index_and_offset(absolute_times, time_zero)
        return time_zero_offset, end_time, time_zero, time_offset, index

NXevent_data_names = ['event_id', 'event_index', 'event_time_offset', 'event_time_zero']
NXevent_data_dtypes = ['i4', 'i8', 'f4', 'f8']

def make_bank(parent, chunk=None, compression=None):
    event_data = parent.create_group("event_data")
    event_data.attrs.create('NX_class', 'NXevent_data', None, dtype='<S12')
    for name, dtype in zip(NXevent_data_names, NXevent_data_dtypes):
        event_data.create_dataset(name, (0,), dtype=dtype, chunks=chunk, compression=compression, maxshape=(None,))
    event_data['event_time_offset'].attrs['units'] = 'microsecond'
    event_data['event_time_zero'].attrs['units'] = 'second'

    return event_data

def link_bank_into_root(root, target, nx_detector_name):
    root["{}_event_data".format(nx_detector_name)] = target

if __name__ == '__main__':
    copyfile(args.reference_filename, args.output_filename)

    with h5py.File(args.input_filename, 'r') as input_file, h5py.File(args.output_filename, 'a') as output_file:
        first_bank = True

        for src, dst, id_offset in zip(input_config['event_data'], output_config['NXdetector'], output_config['event_id_offset']):
            src_data = input_file[input_config['root']][src]
            dst_root = output_file[output_config['root']]
            dst_NXdetector = dst_root[output_config['NXinstrument']][dst]

            make_bank(dst_NXdetector)
            dst_data = dst_NXdetector['event_data']
            link_bank_into_root(dst_root, dst_data, dst)

            print('Loading event data from group {}, writing into {}'.format(src, dst))

            # event_index is just a monotonically increasing integer, not needed.
            # event_index = data['event_index']

            # event_id is 16 bit x pos and 16 bit y pos -> convert
            event_id = convert_id(src_data['event_id'], id_offset)

            # event_time_offset is zero or amplitude of pulse, not needed.
            # event_time_offset = data['event_time_offset']

            # event_time_zero is the absolute time stamp of each event.
            absolute_times = src_data['event_time_zero']
            time_zero_offset, end_time, time_zero, time_offset, index = convert_time(absolute_times)

            dst_data['event_id'].resize(len(event_id), axis = 0)
            dst_data['event_id'][:] = event_id

            dst_data['event_time_zero'].resize(len(time_zero), axis = 0)
            dst_data['event_time_zero'][:] = time_zero
            dst_data['event_time_zero'].attrs['offset'] = time_zero_offset

            dst_data['event_time_offset'].resize(len(time_offset), axis = 0)
            dst_data['event_time_offset'][:] = time_offset

            dst_data['event_index'].resize(len(index), axis = 0)
            dst_data['event_index'][:] = index

            if first_bank:
                dt = h5py.special_dtype(vlen=str)
                dst_root.create_dataset('start_time', (1,), dtype=dt, maxshape=(None,))
                dst_root['start_time'][0] = time_zero_offset
                dst_root.create_dataset('end_time', (1,), dtype=dt, maxshape=(None,))
                dst_root['end_time'][0] = end_time
                first_bank = False
