import numpy as np
import h5py # for header_info

class CatalogueReader:

    def __init__(self, fname, only_coords=False, has_header=True, delim=","):
        """
            fname = filename of catalogue, csv file format
            only_coords = does the csv file only contain the coordinates or also other
                          information? will choose last 3 columns as coords
            has_header = does the csv file have a header?
            delim = delimiter in csv
        """
        self._fname = fname
        self._only_coords = only_coords
        self._header = has_header
        self._delim = delim
        self._data = None

        # extract location of EAGLE data
        # e.g. /net/astrogate/export/astrodata/EAGLE_snapshots/RefL0025N0376/snap_012_z003p017.0.hdf5
        with open(self._fname, 'r') as f:
            self._dataloc = f.readline().strip('# \n')

    def snap_loc(self):
        return self._dataloc

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
        # is always minimum 1 since we need to skip the commentary line
        if self._header:
            skip = 2 + startrow
        else:
            skip = 1 + startrow

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

        if not self._only_coords:
            # extract last 3 columns
            data = data[:,-3:]

        # save data
        self._data = data

        # return to user
        return data
