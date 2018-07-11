# Overview

These tools allow extraction of particles in the vicinity of filaments, as well as 
the visualisation of their properties, such as density _p_, temperature _T_, mass, metallicity, and so on.
Currently this is available in form of a histogram or a _p - T_ phase diagram.

## Approach

The filament catalogue (see `FILAMENT_CATALOGUE`) provides coordinate samples of the filament structure.
Then, for each sample point, we set up a cube with certain box length and the sample point as centre and extract
all particles within this cube. 
These particles are defined to be part of the filament and it is the properties of these particles that are used 
for visualisation.

## Requirements

Note, that this repository does not contain the particle data. 
It provides the filament catalogue but requires data from a cosmological simulation (e.g. EAGLE).

# Code

## Filament extraction

The class `CatalogueReader` allows to open the filament catalogue and extract the filament skeleton.
This is then fed to an instance of `Filaments` which calculates the region of the filament and fetches the particles
from the cosmological simulation.

## Property visualisation

For visualisation the class `Visualiser` is provided. 
See the file `main.py` for an example on how to use it.

## Usage

The driver is `main.py`, see there on how to use this tool.
