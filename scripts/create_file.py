import argparse
import h5py
import numpy

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-p", "--path", type=str, default='', help='Output path used to prefix filename.')
parser.add_argument("-f", "--filename", type=str, default='', help='Output filename.')
parser.add_argument("-n", "--num-elements", type=int, default=32*1024*1024, help='Number of elements to write.')
parser.add_argument("-c", "--chunk", type=int, default=None, help='Chunk size.')
parser.add_argument("-g", "--gzip", action='store_true', help='Enable GZIP compression.')
args = parser.parse_args()

def get_chunking(args):
    if args.chunk is None:
        return None
    return (args.chunk,)

def get_compression(args):
    if args.gzip:
        return "gzip"
    return None

def get_path(args):
    path = args.path
    if len(path) > 0:
        if not path.endswith('/'):
            return "{}/".format(path)
    return path

def get_filename(args):
    filename = args.filename
    if len(filename) > 0:
        if filename.lower().endswith(('.h5', '.hdf5', '.nxs')):
            return path + filename
        return path + filename + ".hdf5"
    return "{}num-element={}_chunk={}_compress={}.hdf5".format(get_path(args), args.num_elements, args.chunk, get_compression(args))

filename = get_filename(args)
count = args.num_elements
chunk = get_chunking(args)
compression = get_compression(args)

f = h5py.File(filename, "w")
dset = f.create_dataset("entry/bank1_events/event_id", (count,), dtype='i', data=numpy.arange(count), chunks=chunk, compression=compression)
