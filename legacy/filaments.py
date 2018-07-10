import numpy as np
import matplotlib.pyplot as plt
from read_eagle import EagleSnapshot
from read_header import read_header
import h5py

from catalogue_reader import CatalogueReader

"""
  [x] Cleanup
  [x] Use coordinates of filament catalogue instead of galaxies
  [ ] Option to remove galaxies from coordinates in filament catalogue
  [ ] Find good load_region_length
  [x] Extract different properties than densities
  [ ] Compare DM densities to baryonic densities
"""

"""
    Global variables
"""
DATADIR = '/net/astrogate/export/astrodata/EAGLE_snapshots/RefL0025N0376/snapshot_012_z003p017/'
SNAPSHOT = 'snap_012_z003p017.0.hdf5'

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

        # Load information from the header
        self._a, self._h, self._boxsize = catalogue.header_info()
        self._dataloc = catalogue.snap_loc()

        # Load centres from filament catalogue
        self._centres = catalogue.load()
        if not (catalogue_index is None):
            self._centres = self._centres[catalogue_index, :]

            # if size is one we need to wrap a list around it for it to be
            # usable in read_particles
            print self._centres.shape
            if len(self._centres.shape) == 1:
                print "fixing"
                self._centres = self._centres.reshape((1,3))

        print "Centres"
        print self._centres
        print self._centres[0]
        # Load data
        self.read_particles()


    def read_particles(self):
        """
            Load particles from simulation
        """

        data = {}

        # Initialize read_eagle module.
        eagle_data = EagleSnapshot(self._dataloc)

        for centre in self._centres:
            # Put centre into cMpc/h units.
            centre *= self._h
            print centre

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

class Visualiser:
    def __init__(self):
        pass

    def hist(self, filaments, att, title="", labels="Data", saveas="_hist.png"):
        """
            filaments: Instance of class Filaments or list of it
            att: attribute to plot
        """
        plt.figure()
        if isinstance(filaments, list):
            if not isinstance(labels, list):
                print "Labels should be a list too!"
                print "This might give strange results."
                # TODO do a fix

            for l, f in zip(labels, filaments):
                # Attribute check
                if not (att in f.atts()):
                    print f, "does not contain the requested attribute."
                    print "Available attributes:"
                    print f.atts()
                    print "Skipping this filament for plotting."

                plt.hist(np.log10(f.data()[att]), bins=50, label=l, normed=True)

        else:
            if not (att in filaments.atts()):
                print "Filaments object does not contain the requested attribute."
                print "Available attributes:"
                print f.atts()
                plt.close()
                return

            plt.hist(np.log10(filaments.data()[att]), bins=50, label=labels, normed=True)

        plt.legend(loc="best")
        plt.minorticks_on()
        plt.title(title)
        plt.xlabel("log10" + att)
        plt.tight_layout()
        plt.savefig(att + saveas)
        plt.close()

if __name__ == '__main__':
    attributes = ['Density', 'Temperature', 'StarFormationRate']
    catalogue = CatalogueReader("FILAMENT_CATALOGUE/s12.csv")

    targets = ["Baryon"]

    if "Filament" in targets:
        # Filament
        region_length = 1

        fil = Filaments(attributes, catalogue, region_length)

        # Plot
        fil.hist('Density', title='Density histogram, RL = ' + str(region_length))
        fil.hist('Temperature', title='Temperature histogram, RL = ' + str(region_length))

    if "Full" in targets:
        region_length = 25 # full box
        catalogue_index = 0 # need only one location
        fil = Filaments(attributes, catalogue, region_length, catalogue_index)

        # Plot
        fil.hist('Density', title='Density histogram, RL = ' + str(region_length), saveas="_hist_full.png")
        fil.hist('Temperature', title='Temperature histogram, RL = ' + str(region_length), saveas="_hist_full.png")

    if "Combined" in targets:
        # Filament
        fil01 = Filaments(attributes, catalogue, region_length=0.1)
        fil02 = Filaments(attributes, catalogue, region_length=0.2)
        fil05 = Filaments(attributes, catalogue, region_length=0.5)
        fil08 = Filaments(attributes, catalogue, region_length=0.8)
        fil1 = Filaments(attributes, catalogue, region_length=1)
        fil2 = Filaments(attributes, catalogue, region_length=2)
        all = Filaments(attributes, catalogue, region_length=25, catalogue_index=0)

        filaments = [all, fil2, fil1, fil08, fil05, fil02, fil01]
        label = ["Full cube (25 Mpc)", "2 Mpc", "1 Mpc", "0.8 Mpc", "0.5 Mpc", "0.2 Mpc", "0.1 Mpc"]

        screen = Visualiser()
        screen.hist(filaments, 'Density', title="Different cube sizes", labels=label, saveas="_hist_compare.png")
        screen.hist(filaments, 'Temperature', title="Different cube sizes", labels=label, saveas="_hist_compare.png")

    if "Baryon" in targets:
        baryon_ctlg = CatalogueReader("/scratch/jgacon/DisPerSE/EAGLE/BARYONS/REFL0012N0188/FILAMENT/s3_baryons.csv")
        fil = Filaments(attributes, baryon_ctlg, region_length=0.01)

        screen = Visualiser()
        screen.hist(fil, 'Density', title='Baryonic filaments', saveas="_hist_bary.png")
