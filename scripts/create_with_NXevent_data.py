import argparse
import h5py
import numpy

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-P", "--path", type=str, default='', help='Output path used to prefix filename.')
parser.add_argument("-f", "--filename", type=str, default='', help='Output filename.')
parser.add_argument("-e", "--events", type=float, default=1e6, help='Number of events to write (approximate).')
parser.add_argument("-c", "--chunk", type=int, default=256*1024, help='Chunk size.')
parser.add_argument("-g", "--gzip", action='store_true', help='Enable GZIP compression.')
parser.add_argument("-b", "--banks", type=int, default=7, help='Number of banks.')
parser.add_argument("-p", "--pixels", type=int, default=10000, help='Number of detector pixels (per bank).')
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
    return "{}events-{}_banks-{}_pixels-{}_chunk-{}_compress-{}.hdf5".format(get_path(args), int(args.events), args.banks, args.pixels, args.chunk, get_compression(args))

NXevent_data_names = ['event_index', 'event_time_zero', 'event_id', 'event_time_offset']
NXevent_data_dtypes = ['i4', 'i8', 'i4', 'i4']

def make_bank(instrument, index, chunk, compression):
    event_data = instrument.create_group("events-{}".format(index))
    event_data.attrs['NXclass'] = 'NXevent_data'
    for name, dtype in zip(NXevent_data_names, NXevent_data_dtypes):
        event_data.create_dataset(name, (0,), dtype=dtype, chunks=chunk, compression=compression, maxshape=(None,))
    return event_data

def extend_bank(bank, delta):
    NXevent_data_increments = [1, 1, delta, delta]
    for name, delta in zip(NXevent_data_names, NXevent_data_increments):
        bank[name].resize(bank[name].size + delta, axis = 0)

filename = get_filename(args)

f = h5py.File(filename, "w")
entry = f.create_group("entry")
entry.attrs['NXclass'] = 'NXentry'
instrument = entry.create_group("instrument")
instrument.attrs['NXclass'] = 'NXinstrument'

chunk = get_chunking(args)
compression = get_compression(args)

for bank in range(args.banks):
    event_data = make_bank(instrument, bank, chunk, compression)

numpy.random.seed(0)
bank_load_factor = numpy.random.exponential(size=args.banks)
bank_size = args.pixels # per bank
bank_offset = 10 * bank_size # spread factor to make event_id -> workspace index translation less trivial (but all pixels within bank are assumed to be contiguous for now)

max_events = int(args.events)
total_events = 0
pulse = 0
while total_events < max_events:
    n_event_base_scale = 10000 # rough number of events going into each bank per pulse
    n_event_base = n_event_base_scale*abs(numpy.random.normal(loc=1.0))
    for i, bank in enumerate(instrument.values()):
        n_event = int(bank_load_factor[i] * n_event_base)
        total_events += n_event
        extend_bank(bank, n_event)
        bank['event_index'][-1:] = bank['event_id'].size
        bank['event_time_zero'][-1:] = pulse * 100000
        bank['event_id'][-n_event:] = numpy.random.randint(i*bank_offset, i*bank_offset+bank_size, size=n_event, dtype=numpy.int32)
        bank['event_time_offset'][-n_event:] = numpy.random.randint(10000, 20000, size=n_event, dtype=numpy.int32)
    pulse += 1
