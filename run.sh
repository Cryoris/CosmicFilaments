#! /bin/bash

# how many filament centres are being processed at once
stepsize=50000
# Galaxy catalogue
#filament_catalogue="FILAMENT_CATALOGUE/s12.csv"
# Baryon catalogue
filament_catalogue="/scratch/jgacon/DisPerSE/EAGLE/BARYONS/REFL0012N0188/FILAMENT/s3_baryons.csv"
ncentres=$(wc -l < $filament_catalogue)
ncentres=$(expr $ncentres - 2)

echo "Number of centres: $ncentres"
echo "Stepsize: $stepsize"

count=0
while [ $count -lt $ncentres ]
do
  idx_low=$count
  count=$(expr $count + $stepsize)
  idx_high=$count
  # actual command that runs
  # TODO job system
  python main.py $idx_low $idx_high
done
