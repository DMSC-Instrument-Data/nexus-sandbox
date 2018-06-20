# Take an existing Nexus event file and increase the number of events in it

import argparse
import h5py
import numpy

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-f", "--filename", type=str, default='', help='File to grow. Modified in place.')
parser.add_argument("-s", "--scale", type=float, default=2, help='Scale factor to multiply the number of events.')
parser.add_argument("-r", "--root", type=str, default='entry', help='Name of the root entry containing the NXevent_data entries.')
args = parser.parse_args()

def grow(filename):
    with h5py.File(filename, 'a') as f:
        scale = args.scale
        if scale < 1:
            print('Warning: Reducing event count, but chunking is left untouched. If the chunk size is large the resulting file may actually not shrink. If that is the case, use something like `h5repack -i {0}.nxs -o {0}-rechunked.nxs -l CHUNK=N` to produce a smaller file.'.format(filename[:-4]))
        for bank, event_data in f[args.root].items():
            if not bank.endswith('_events'):
                continue
            print("Found bank {}".format(bank))
            if scale > 1:
                event_data['event_index'][:] = event_data['event_index'][:] * scale
            else:
                pulses = event_data['event_index'].size
                event_data['event_index'].resize(scale*pulses, axis = 0)
                event_data['event_time_zero'].resize(scale*pulses, axis = 0)
            size = event_data['event_id'].size
            event_data['event_id'].resize(scale*size, axis = 0)
            event_data['event_time_offset'].resize(scale*size, axis = 0)
            for s in range(int(scale)):
                event_data['event_id'][s*size:(s+1)*size] = event_data['event_id'][:size]
                event_data['event_time_offset'][s*size:(s+1)*size] = event_data['event_time_offset'][:size]

grow(args.filename)
