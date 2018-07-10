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


    def load(self):
        """
            read data and save to self._data and return it
        """
        # is always minimum 1 since we need to skip the commentary line
        if self._header:
            skip = 2
        else:
            skip = 1

        print "Loading data from", self._fname,"..."
        data = np.loadtxt(self._fname, delimiter=self._delim, skiprows=skip)
        print "Loaded", data.shape[0], "particles."

        if not self._only_coords:
            # extract last 3 columns
            data = data[:,-3:]

        # save data
        self._data = data

        # return to user
        return data
