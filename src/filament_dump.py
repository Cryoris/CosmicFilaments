import numpy as np
from shutil import copyfile
from subprocess import Popen

class FilamentDump:
    def __init__(self):
        pass

    def dump(self, fdict, outfile):
        """
            Dump filament dictionary to csv
            -> keys as rownames and data as columns
        """
        header = str(fdict.keys()).strip("[]").replace("'","")

        n_att = len(fdict) # Number of attributes / columns
        n_points = fdict.values()[0].size # Number of points

        data = np.empty((n_points, n_att))

        for i, vals in enumerate(fdict.values()):
            data[:,i] = vals

        np.savetxt(outfile, data, header=header, delimiter=",")

    def gather(self, files, outfile):
        # copy first, then append data from others
        if len(files) > 0:
            copyfile(files[0], outfile)

        if not len(files) > 1:
            return

        for f in files[1:]:
            # Skip first row (= rownames and append rest)
            bash = "tail -n +2 {} >> {}".format(f, outfile)
            Popen(bash, shell=True)
