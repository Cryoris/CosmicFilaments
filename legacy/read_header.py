import h5py

def read_header():
    """ Read various attributes from the header group. """
    f       = h5py.File('/net/astrogate/export/astrodata/EAGLE_snapshots/RefL0025N0376/snapshot_012_z003p017/snap_012_z003p017.0.hdf5', 'r')
    a       = f['Header'].attrs.get('Time')         # Scale factor.
    h       = f['Header'].attrs.get('HubbleParam')  # h.
    boxsize = f['Header'].attrs.get('BoxSize')      # L [Mph/h].
    f.close()

    return a, h, boxsize
