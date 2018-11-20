import argparse
import h5py
import numpy

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i", "--input-filename", type=str, help='Input file to convert.')
parser.add_argument("-r", "--reference-filename", type=str, help='Input file used as a reference for the file structure.')
#parser.add_argument("-s", "--scale", type=float, default=2, help='.')
#parser.add_argument("-r", "--root", type=str, default='entry', help='Name of the root entry containing the NXevent_data entries.')
args = parser.parse_args()

input_config = dict(root='entry-01', event_data=['Delayline_events'], Nx=8, Ny=16)

def position_to_index(pos, count):
    uint_max = 2^32-1
    # What if count (Nx or Ny) does not divide uint_max?
    return np.floor_divide(pos/(uint_max//count))


def convert_id(event_id):
    x = (event_id >> 16) & 0xffff
    y = event_id & 0xffff
    Nx = input_config['Nx']
    Ny = input_config['Ny']
    return position_to_index(x, Nx) + Nx * position_to_index(y, Ny)

if __name__ == '__main__':
    with h5py.File(args.input_filename, 'r') as input_file, h5py.File(args.reference_filename, 'r') as reference:
        for group, data in input_file[input_config['root']].items():
            if group in input_config['event_data']:
                print('Loading event data from group {}'.format(group))

                event_index = data['event_index']
                event_id = data['event_id']
                event_offset = data['event_offset']
                event_time_zero = data['event_time_zero']


