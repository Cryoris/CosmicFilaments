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

class FilamentProperties:
    """
        Class to handle read-in and visualisation of particle properties near filaments
    """
    def __init__(self, gn, sgn, centres, load_region_length=1):
        """
            gn: GroupNumber
            sgn: SubGroupNumber
            centre: centre of region of particles that we look at
            load_region_length: cube length of particle region
        """

        # Load information from the header.
        self.a, self.h, self.boxsize = read_header()

        # Load data.
        self.gas = self.read_particles(0, gn, sgn, centres, load_region_length)

        # Plot.
        #self.tempdensity()
        self.hist('Density', title='Density histogram, RL = ' + str(load_region_length))
        self.hist('Temperature', title='Temperature histogram, RL= ' + str(load_region_length))
        #self.hist('StarFormationRate')

    def read_particles(self, itype, gn, sgn, centres, load_region_length):
        """
            For a given galaxy (defined by its GroupNumber and SubGroupNumber)
            extract the Temperature, Density and StarFormationRate of all gas particles
            using the read_eagle routine. Conversion factors are still loaded directly
            from the hdf5 files.

            itype: Particle Type (0 for PartType0 = Baryons)
            gn: GroupNumber
            sgn: SubGroupNumber
            centres: centres of region of particles that we look at
                     -> array of dim: ncentres x 3
            load_region_length: cube length of particle region
                     -> constant number
        """

        data = {}

        # Initialize read_eagle module.
        eagle_data = EagleSnapshot(DATADIR + SNAPSHOT)

        for centre in centres:
            # Put centre into cMpc/h units.
            print centre
            centre *= self.h

            # Select region to load, a 'load_region_length' cMpc/h cube centred on 'centre'.
            region = np.array([
                (centre[0]-0.5*load_region_length), (centre[0]+0.5*load_region_length),
                (centre[1]-0.5*load_region_length), (centre[1]+0.5*load_region_length),
                (centre[2]-0.5*load_region_length), (centre[2]+0.5*load_region_length)
            ])
            eagle_data.select_region(*region)

            # Load data using read_eagle, load conversion factors manually.
            f = h5py.File(DATADIR + SNAPSHOT, 'r')
            for att in ['Temperature', 'Density', 'StarFormationRate']:
                tmp  = eagle_data.read_dataset(itype, att)
                cgs  = f['PartType%i/%s'%(itype, att)].attrs.get('CGSConversionFactor')
                aexp = f['PartType%i/%s'%(itype, att)].attrs.get('aexp-scale-exponent')
                hexp = f['PartType%i/%s'%(itype, att)].attrs.get('h-scale-exponent')
                data[att] = np.multiply(tmp, cgs * self.a**aexp * self.h**hexp, dtype='f8')
            f.close()

            print type(data)
            print type(data[data.keys()[0]])
            print "Size of data:", len(data[data.keys()[0]])

            # Mask to selected GroupNumber and SubGroupNumber.
            """
            mask = np.logical_and(data['GroupNumber'] == gn, data['SubGroupNumber'] == sgn)
            for att in data.keys():
                data[att] = data[att][mask]
            """

        return data



    def read_galaxy(self, itype, gn, sgn, centre, load_region_length):
        """
            For a given galaxy (defined by its GroupNumber and SubGroupNumber)
            extract the Temperature, Density and StarFormationRate of all gas particles
            using the read_eagle routine. Conversion factors are still loaded directly
            from the hdf5 files.

            itype: Particle Type (0 for PartType0 = Baryons)
            gn: GroupNumber
            sgn: SubGroupNumber
            centre: centre of region of particles that we look at
            load_region_length: cube length of particle region
        """

        data = {}

        # Initialize read_eagle module.
        eagle_data = EagleSnapshot(DATADIR + SNAPSHOT)

        # Put centre into cMpc/h units.
        centre *= self.h

        # Select region to load, a 'load_region_length' cMpc/h cube centred on 'centre'.
        region = np.array([
            (centre[0]-0.5*load_region_length), (centre[0]+0.5*load_region_length),
            (centre[1]-0.5*load_region_length), (centre[1]+0.5*load_region_length),
            (centre[2]-0.5*load_region_length), (centre[2]+0.5*load_region_length)
        ])
        eagle_data.select_region(*region)

        # Load data using read_eagle, load conversion factors manually.
        f = h5py.File(DATADIR + SNAPSHOT, 'r')
        for att in ['GroupNumber', 'SubGroupNumber', 'Temperature', 'Density', 'StarFormationRate']:
            tmp  = eagle_data.read_dataset(itype, att)
            cgs  = f['PartType%i/%s'%(itype, att)].attrs.get('CGSConversionFactor')
            aexp = f['PartType%i/%s'%(itype, att)].attrs.get('aexp-scale-exponent')
            hexp = f['PartType%i/%s'%(itype, att)].attrs.get('h-scale-exponent')
            data[att] = np.multiply(tmp, cgs * self.a**aexp * self.h**hexp, dtype='f8')
        f.close()

        print "Size of data:", len(data)

        # Mask to selected GroupNumber and SubGroupNumber.
        """
        mask = np.logical_and(data['GroupNumber'] == gn, data['SubGroupNumber'] == sgn)
        for att in data.keys():
            data[att] = data[att][mask]
        """
        return data

    def hist(self, att='Density', title='', saveas='_hist.png'):
        """
            Plot histogram
        """
        plt.figure()

        plt.hist(np.log10(self.gas[att]), bins=20)

        # Save plot.
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
        plt.figure()

        # Plot currently star forming gas red.
        mask = np.where(self.gas['StarFormationRate'] > 0)
        plt.scatter(np.log10(self.gas['Density'][mask]), np.log10(self.gas['Temperature'][mask]),
            c='red', s=3, edgecolor='none')

        # Plot currently non star forming gas blue.
        mask = np.where(self.gas['StarFormationRate'] == 0)
        plt.scatter(np.log10(self.gas['Density'][mask]), np.log10(self.gas['Temperature'][mask]),
            c='blue', s=3, edgecolor='none')

        # Save plot.
        plt.minorticks_on()
        plt.ylabel('log10 Temperature [K]'); plt.xlabel('log10 Density [g/cm**3]')
        plt.tight_layout()
        plt.savefig('PhaseDiagram_s12.png')
        plt.show()
        plt.close()

if __name__ == '__main__':
    # Centre is the COP for GN=1 SGN=0 taken from the database.
    catalogue = CatalogueReader("FILAMENT_CATALOGUE/s12.csv")
    centres = catalogue.get()
    ngals = -1 # -1 == all galaxies
    if ngals > -1:
        centres = centres[:ngals,:]


    #centre = np.array([12.08808994,4.47437191,1.41333473])  # cMpc
    x = FilamentProperties(1, 0, centres, 0.5)
