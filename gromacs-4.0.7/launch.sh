#!/bin/bash

# An example to launch a simulation

# Change the path to correspond to your path to force field files
source /usr/local/gromacs-4.0.7/bin/GMXRC
export GMXLIB=/MY_PATH/CiEj_2016H66/gromacs-4.0.7/ffG2016H66
export GMXDATA=/MY_PATH/CiEj_2016H66/gromacs-4.0.7/ffG2016H66

# .itp file must be in the same repertory as topol.top file
grompp -f mdin.mdp -c start.gro -p topol.top -o md
mdrun -v -deffnm md

