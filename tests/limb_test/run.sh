#! /bin/bash

# Set environment...
export LD_LIBRARY_PATH=../../libs/build/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=4
export LANG=C
export LC_ALL=C

# Setup...
jurassic=../../src

# Create directory...
rm -rf data && mkdir -p data || exit

# Create atmospheric data file...
$jurassic/climatology limb.ctl data/atm.tab

# Create observation geomtry...
$jurassic/limb limb.ctl data/obs.tab

# Call forward model...
$jurassic/formod limb.ctl data/obs.tab data/atm.tab data/rad.tab TASK time

# Test CGA...
$jurassic/formod limb.ctl data/obs.tab data/atm.tab data/rad_cga.tab FORMOD 0

# Test FOV...
$jurassic/formod limb.ctl data/obs.tab data/atm.tab data/rad_fov.tab FOV fov.tab

# Compute kernel...
$jurassic/kernel limb.ctl data/obs.tab data/atm.tab data/kernel.tab

# Compare files...
echo -e "\nCompare results..."
error=0
for f in $(ls data.ref/*.tab) ; do
    diff -q -s data/"$(basename "$f")" "$f" || error=1
done
exit $error
