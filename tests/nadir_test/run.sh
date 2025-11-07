#! /bin/bash

# Set environment...
export LD_LIBRARY_PATH=../../libs/build/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=4
export LANG=C
export LC_ALL=C

# Setup...
jurassic=../../src

# Create directory...
rm -rf data && mkdir -p data

# Create atmospheric data file...
$jurassic/climatology nadir.ctl data/atm.tab

# Create observation geomtry...
$jurassic/nadir nadir.ctl data/obs.tab

# Call forward model...
$jurassic/formod nadir.ctl data/obs.tab data/atm.tab data/rad.tab TASK time

# CGA test...
$jurassic/formod nadir.ctl data/obs.tab data/atm.tab data/rad_cga.tab FORMOD 0

# Compute kernel...
$jurassic/kernel nadir.ctl data/obs.tab data/atm.tab data/kernel.tab

# Compare files...
echo -e "\nCompare results..."
error=0
for f in $(ls data.ref/*.tab) ; do
    diff -q -s data/"$(basename "$f")" "$f" || error=1
done
exit $error
