# Adapted from the original by Matthew D Jones to be found at:
# https://github.com/ScreamingUdder/isis_nexus_streamer_for_mantid/blob/master/data/nexus_recompress_all.py
# Changes should make the output equivalent to the original NeXus files such that Mantid can read it. Changes include:
# - fixed data types
# - copy attributes
# - set up links
# TODO: Cleanup

# export HDF5_PLUGIN_PATH=/home/simon/software/hdf5-blosc/build

import h5py
import tables
import os

"""
Read data from nexus file and write it to a new BLOSC compressed one.
"""

original_file = 'PG3_4871_event.nxs'

def h5py_dataset_iterator(g, prefix=''):
    for key in g.keys():
        item = g[key]
        path = '{}/{}'.format(prefix, key)
        if isinstance(item, h5py.Dataset): # test for dataset
            yield (path, item)
        elif isinstance(item, h5py.Group): # test for group (go down)
            yield (path, item)
            yield from h5py_dataset_iterator(item, path)

def write_file(output_file, compress_type=32001, compress_opts=None):
    with h5py.File(output_file, 'w') as f_write:
        with h5py.File(original_file, 'r') as f_read:

            linked = []

            def write_in_new_file(name):
                dset = f_read[name]
                if 'target' in dset.attrs and name != dset.attrs['target']:
                    linked.append((str(dset.attrs['target'],'utf-8'),name))
                if isinstance(dset, h5py.Dataset):
                    try:
                        output_dataset = f_write.create_dataset(name, dset[...].shape, dtype=dset.dtype, compression=compress_type,
                                                                compression_opts=compress_opts)
                        output_dataset[...] = dset[...]
                    except TypeError:
                        # probably this is a scalar dataset so just write it without compression
                        print(name)
                        f_write[name] = dset[...]
                else:
                    f_write.create_group(name)
                # Copy all attributes
                for key,value in dset.attrs.items():
                    f_write[name].attrs[key] = value

            f_read.visit(write_in_new_file)

            for path, target in linked:
                if not path in f_write:
                    try:
                        f_write[path] = f_write[target]
                        print('Set up link from {} to {}'.format(path, str(target)))
                    except KeyError:
                        print('Cannot find target for {}'.format(path))

            for (path, dset) in h5py_dataset_iterator(f_read):
                if not path in f_write:
                    #f_write['{}/time'.format(name)] = f_write['entry/DASlogs/frequency/time']
                    try:
                        f_write[path] = f_write[str(dset.attrs['target'],'utf-8')]
                        print('Set up link from {} to {}'.format(path, str(dset.attrs['target'])))
                    except KeyError:
                        print('Cannot find target for {}'.format(path))


if __name__ == '__main__':
    gzip_file_name = 'gzip.nxs'
    blosc_file_name = 'blosc.nxs'

    #write_file(gzip_file_name, compress_type='gzip')
    write_file(blosc_file_name, compress_type=32001)  # 32001 is the Blosc filter
