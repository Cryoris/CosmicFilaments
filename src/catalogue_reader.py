import numpy as np
import h5py # for header_info

"""
    About the structure of the filament catalogues:
    The number of commentary lines (4) and the header is mandatory!

    # /location/of/EAGLE/snapshot/that/this/catalogue/uses/snapshot.hdf5
    # cube size [Mpc]
    # number of particles/resolution [integer]
    # particle type used & other metadata [string]
    colname_1, colname_2, ..., colname_n-3, x, y, z
    DATA_1, DATA_2, ..., DATA_n-3, x coord, y coord, z coord
    DATA_1, DATA_2, ..., DATA_n-3, x coord, y coord, z coord
    DATA_1, DATA_2, ..., DATA_n-3, x coord, y coord, z coord
    etc.

    We don't care about the first n-3 columns, they might contain additional
    data, we only read out the last 3 columns. Of course the catalogue can
    only contain 3 columns (the coordinates), but in case more data is
    supposed to be saved this is possible.
"""

class CatalogueReader:

    def __init__(self, fname, delim=","):
        """
            fname = filename of catalogue, csv file format
            has_header = does the csv file have a header?
            delim = delimiter in csv
        """
        self._fname = fname
        self._delim = delim
        self._data = None

        # extract location of EAGLE data
        # e.g. /net/astrogate/export/astrodata/EAGLE_snapshots/RefL0025N0376/snap_012_z003p017.0.hdf5
        # and cube size in Mpc
        with open(self._fname, 'r') as f:
            self._dataloc = f.readline().strip('# \n')
            self._cube_length = float(f.readline().strip('# \n'))
            self._num_particles = int(f.readline().strip('# \n'))
            self._metadata = f.readline().strip('# \n')

    def snap_loc(self):
        return self._dataloc

    def cube_length(self):
        return self._cube_length

    def num_particles(self):
        """
            Note that this is not the number of centres, but the number of
            particles used in the creation of the filament!
            This is a measure of how precise/coarse grained the catalogue is.
        """
        return self._num_particles

    def metadata(self):
        """
            Content of 4th commentary line.
            Usually what kind of particles have been used for creation of the catalogue
             (--> Baryons, Galaxies, Dark Matter)
        """
        return self._metadata

    def num_centres(self):
        # minus 2 since first line is location of the snapshot and
        # second is the rownames
        return sum(1 for line in open(self._fname)) - 2

    def header_info(self):
        """
            Read various attributes from the header group.
        """
        f = h5py.File(self._dataloc, 'r')
        a = f['Header'].attrs.get('Time')         # Scale factor.
        h = f['Header'].attrs.get('HubbleParam')  # h.
        boxsize = f['Header'].attrs.get('BoxSize')      # L [Mph/h].
        f.close()

        return a, h, boxsize

    def load(self, startrow=0, endrow=None):
        """
            read data and save to self._data and return it
            startrow: in which row to start reading the data (first row = 0)
            endrow: last row (0-based count)
        """
        # number of lines to skip: 4 (commentary) + 1 (header) = 5
        skip = 5 + startrow

        print "Loading data from", self._fname,"..."
        if endrow is None:
            data = np.genfromtxt(self._fname, delimiter=self._delim, skip_header=skip)
        else:
            maxrows = endrow - startrow
            data = np.genfromtxt(self._fname, delimiter=self._delim, skip_header=skip, max_rows=maxrows)
        print "Loaded", data.shape[0], "particles."

        # ensure 2d (matrix form)
        if len(data.shape) == 1:
            data = data.reshape((1, data.size))

        # extract last 3 columns (also valid if only 3 exist)
        # and save data
        self._data = data[:,-3:]

        # return to user
        return self._data
