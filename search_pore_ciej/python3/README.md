This directory contains the Python 3 version of `search_pore_ciej.py` with added features.

-------------------------
search_pore_com_ciej.py is a program to find pores in a bilayer and compute their center of mass (COM).  
Author: Caroline Senac (caroline.senac@upmc.fr) and Maya Zygadlo (maya.zygadlo@sorbonne-universite.fr)  
Version: V4 (07/2024)  

This program cuts the bilayer in a grid of 200*200 cells (can be modified by the user)
and searches for empty cells. The connected empty cells are grouped in different pores. 
Then, their center of mass is computed taking into account the periodic boundary conditions.


Prerequisites:
- Python 3
- Python libraries : scipy, numpy, mdtraj

Usage:  
```
python search_pore_com_ciej.py -f traj.xtc -c topol.gro -n index.ndx -p list_atoms.txt [-col 200] [-row 200] [-step 500000]
```

Mandatory arguments:  
  -f TRAJECTORY, --trajectory TRAJECTORY    <.xtc/.gro> trajectory file  
  -c TOPOLOGY, --topology TOPOLOGY    <.gro> topology file  
  -n INDEX, --index INDEX    <.ndx> index file  
  -p LIST_ATOMS, --list_atoms LIST_ATOMS    <.txt> list of atoms  
Optionnal arguments:  
  -h, --help            shows this help message and exit  
  -col COLUMN, --column COLUMN    <int> number of column in grid  
  -row ROW, --row ROW    <int> number of row in grid  
  -step STEP, --step STEP    <int> number of steps between each framein the .xtc file  


Input files:
- A GROMACS trajectory file (traj.xtc). You need to do a "gmx trjconv -pbc atom" on your trajectory to use this program.
- A GROMACS topology file (topol.gro).
- A GROMACS index file (index.ndx). This index must have only one group, which corresponds to the study bilayer (example: aliphatic tails for CiEj).
- A file text of atoms of interest. This file need to be organised in 3 columns separated by spaces.
                                    First column: residue name, second column: atom name, third column type name.
                                    These names must be like in itp file (See the example file: itp_for_pore.txt).

Output:
- pdb files in a repertory named pdb. Each pdb corresponds to a frame and each point in the file is an empty cell.
- pores_in_time.dat: a file with the first and last frame of each pore.
                     Each line goes as follow: pore name;first frame;last frame;number of frame;percentage of frame with this pore.
- sizes.dat: matrix of pore size in nm^2 for each frame.
             Lines: pore (first line is bilayer area in nm^2)
             Columns: frame (first column is pore name)
- pores_size.txt: for each frame the number of pores and their sizes.
- center_of_mass_pore.pdb: a multipdb containing the COM positions of each pore in each frame.


!!! WARNING !!! Delete old output if you rerun the program in the same repertory.


--------------------------------------------------
DISCLAIMER
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
