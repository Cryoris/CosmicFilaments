import sys
from src.filament_extractor import *
from src.visualiser import *

attributes = ['Density', 'Temperature', 'StarFormationRate']
catalogue = CatalogueReader("FILAMENT_CATALOGUE/s12.csv")

# If indices for catalogue data is given define it!
# Use to parallelise data loading
argc = len(sys.argv)
if argc > 3:
   catalogue_index = np.arange(sys.argv[1], sys.argv[2])
elif argc > 2:
    catalogue_index = sys.argv[1]

targets = ["Baryon"]

if "Test" in targets:
    region_length = 0.01 # small test region
    catalogue_index = 0 # need only one location
    fil = Filaments(attributes, catalogue, region_length, catalogue_index)

    # Dump data
    files = ["test{}.dump".format(i) for i in range(1,5)]
    for f in files:
        fil.dump(f)

    fil.gather(files, "test.dump")

    # Plot
    fil.hist('Density', title='Density histogram, RL = ' + str(region_length), saveas="_test.png")
    fil.hist('Temperature', title='Temperature histogram, RL = ' + str(region_length), saveas="_test.png")

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
    fil = Filaments(attributes, baryon_ctlg, region_length=0.01, catalogue_index=catalogue_index)

    fil.dump("fil.{}.dump".format(catalogue_index[0]))

    #screen = Visualiser()
    #screen.hist(fil, 'Density', title='Baryonic filaments', saveas="_hist_bary.png")
