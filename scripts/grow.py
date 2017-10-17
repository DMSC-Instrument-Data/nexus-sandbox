# Take an existing Nexus event file and increase the number of events in it

import argparse
import h5py
import numpy

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-f", "--filename", type=str, default='', help='File to grow. Modified in place.')
parser.add_argument("-s", "--scale", type=int, default=2, help='Scale factor to multiply the number of events.')
args = parser.parse_args()

def grow(filename):
    with h5py.File(filename, 'a') as f:
        for bank, event_data in f['entry'].items():
            if not bank.endswith('_events'):
                continue
            print("Found bank {}".format(bank))
            scale = args.scale
            event_data['event_index'][:] = event_data['event_index'][:] * scale
            size = event_data['event_id'].size
            event_data['event_id'].resize(scale*size, axis = 0)
            event_data['event_time_offset'].resize(scale*size, axis = 0)
            for s in range(scale):
                event_data['event_id'][s*size:(s+1)*size] = event_data['event_id'][:size]
                event_data['event_time_offset'][s*size:(s+1)*size] = event_data['event_time_offset'][:size]

grow(args.filename)
