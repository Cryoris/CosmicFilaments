# General imports
import numpy as np
import h5py
import matplotlib.pyplot as plt

# EAGLE imports
from read_eagle import EagleSnapshot
#from read_header import read_header

# Filament catalogue imports
from catalogue_reader import CatalogueReader
from filament_dump import FilamentDump

class Filaments:
    """
        Class to handle read-in and visualisation of particle properties near filaments
    """
    def __init__(self, attributes, catalogue, region_length=1, catalogue_index=None, part_type=0):
        """
            attributes: which attributes to load?
                        -> e.g. ['Density', 'Temperature']
            catalogue: CatalogueReader with filament catalogue
                       --> see file catalogue_reader.py
            catalogue_index: all particles in filament catalogue (None)
                             or just certain ones (e.g. 1:8)?
            region_length: cube length of particle region
        """

        self._att = attributes
        self._region_length = region_length
        self._itype = part_type
        self._data = None # will be set by read_particles

        # Set instance of FilamentDump, might need it later
        self._dumper = FilamentDump()

        # Load information from the header
        self._a, self._h, self._boxsize = catalogue.header_info()
        self._dataloc = catalogue.snap_loc()

        # Load centres from filament catalogue
        self._centres = catalogue.load()
        if not (catalogue_index is None):
            self._centres = self._centres[catalogue_index, :]

            # if size is one we need to wrap a list around it for it to be
            # usable in read_particles
            if len(self._centres.shape) == 1:
                self._centres = self._centres.reshape((1,3))

        #print "Centres"
        #print self._centres
        # Load data
        self.read_particles()


    def read_particles(self):
        """
            Load particles from simulation
        """

        data = {}

        # Initialize read_eagle module.
        eagle_data = EagleSnapshot(self._dataloc)

        for i, centre in enumerate(self._centres):
            # Put centre into cMpc/h units.
            centre *= self._h
            if i % 10 == 0:
                print "Current box centre (No. {}): is".format(i+1), centre

            # Select region to load, a 'load_region_length' cMpc/h cube centred on 'centre'.
            region = np.array([
                (centre[0]-0.5*self._region_length), (centre[0]+0.5*self._region_length),
                (centre[1]-0.5*self._region_length), (centre[1]+0.5*self._region_length),
                (centre[2]-0.5*self._region_length), (centre[2]+0.5*self._region_length)
            ])
            eagle_data.select_region(*region)

            # Load data using read_eagle, load conversion factors manually.
            f = h5py.File(self._dataloc, 'r')
            for att in self._att:
                tmp  = eagle_data.read_dataset(self._itype, att)
                cgs  = f['PartType%i/%s'%(self._itype, att)].attrs.get('CGSConversionFactor')
                aexp = f['PartType%i/%s'%(self._itype, att)].attrs.get('aexp-scale-exponent')
                hexp = f['PartType%i/%s'%(self._itype, att)].attrs.get('h-scale-exponent')
                data[att] = np.multiply(tmp, cgs * self._a**aexp * self._h**hexp, dtype='f8')
            f.close()

            # Mask to selected GroupNumber and SubGroupNumber.
            """
            mask = np.logical_and(data['GroupNumber'] == gn, data['SubGroupNumber'] == sgn)
            for att in data.keys():
                data[att] = data[att][mask]
            """


        self._data = data

    def atts(self):
        return self._att

    def data(self):
        return self._data

    def dump(self, outfile):
        self._dumper.dump(self._data, outfile)

    def gather(self, files, outfile):
        self._dumper.gather(files, outfile)

    def hist(self, att='Density', title='', saveas='_hist.png'):
        """
            Plot histogram
        """

        if not (att in self._att):
            print "Invalid attribute, didn't load this from the data."
            print "Available attributes:"
            print self._att
            return

        plt.figure()
        plt.hist(np.log10(self._data[att]), bins=20)
        plt.minorticks_on()
        plt.title(title)
        plt.xlabel("log10" + att)
        plt.tight_layout()
        plt.savefig(att + saveas)
        plt.close()

    def tempdensity(self):
        """
            Plot Temperature--Density relation.
        """
        if not ('StarFormationRate' in self._att
                and 'Density' in self._att
                and 'Temperature' in self._att):
            print "Must load 'StarFormationRate', 'Density' and 'Temperature' attribute from data for this function."
            print "Loaded attributes:"
            print self._att
            return

        plt.figure()

        # Plot currently star forming gas red.
        mask = np.where(self._data['StarFormationRate'] > 0)
        plt.scatter(np.log10(self._data['Density'][mask]), np.log10(self._data['Temperature'][mask]),
            c='red', s=3, edgecolor='none')

        # Plot currently non star forming gas blue.
        mask = np.where(self.gas['StarFormationRate'] == 0)
        plt.scatter(np.log10(self._data['Density'][mask]), np.log10(self._data['Temperature'][mask]),
            c='blue', s=3, edgecolor='none')

        # Save plot.
        plt.minorticks_on()
        plt.ylabel('log10 Temperature [K]'); plt.xlabel('log10 Density [g/cm**3]')
        plt.tight_layout()
        plt.savefig('PhaseDiagram_s12.png')
        plt.show()
        plt.close()
