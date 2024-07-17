#!/usr/bin/python
# -*- coding: utf-8 -*-

__author__ = "Caroline SENAC and Maya ZYGADLO"
__date__ = "07/2024"
__version__ = "V4"

import re
import argparse
import sys
import os
import errno

from scipy.ndimage import label, sum
import numpy as np
import mdtraj as md

usageprog = \
"""
    python search_pore_ciej.py -f traj.xtc -c topol.gro -n index.ndx -p list_atoms.txt [-col 200] [-row 200]
"""


def recup_args():
    """
    Get the arguments
    """
    parser = argparse.ArgumentParser(description="Find pore in CiEj bilayer", usage=usageprog)
    parser.add_argument("-f", "--trajectory", help="<.xtc/.gro> trajectory file")
    parser.add_argument("-c", "--topology", help="<.gro> topology file")
    parser.add_argument("-n", "--index", help="<.ndx> index file")
    parser.add_argument("-p", "--list_atoms", help="<.txt> list of atoms")
    parser.add_argument("-col", "--column", type=float, default=200.,\
                        help="<int> number of column in grid. Default : 200")
    parser.add_argument("-row", "--row", type=float, default=200.,\
                        help="<int> number of row in grid. Default : 200")
    parser.add_argument("-step", "--step", type=int, default=500000,\
                        help="<int> number of steps between each framein the .xtc file. Default : 500000")
    args = parser.parse_args()
    return args


def make_dir(dirpath):
    """
    Create a directory
    
    -----------------------------------
    INPUT : 
    dirpath : string
        Path and name of the directory

    """
    try:
        os.makedirs(dirpath)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
        else:
            print(f"Directory {dirpath} already exists!")


def read_ndx(filename):
    """
    Read index file

    -----------------------------------
    INPUT : 
    filename : string
        Name of the index file

    -----------------------------------
    OUTPUT : 
    list
        List of the sorted index of the atom in the entry file
    """
    l_index = []
    regex = re.compile("\[")
    # Read the index file and get the lines
    with open(filename, 'r') as f:
        lines = f.readlines()
    # Get the index of the atoms that are in the index file
    for l in lines:
        # Skip the index names
        if not regex.match(l):
            l_index += [int(i) for i in l.split()]
    # sort the atom index and substract 1 to each
    l_index = np.array(sorted(l_index)) -1
    return l_index


def get_atoms(filename):
    """
    Read atom list file

    -----------------------------------
    INPUT : 
    filename : string
        Name of the atom list file

    -----------------------------------
    OUTPUT : 
    dictionnary
        Dictionnary containing the molecule and general name of a group
        example : {('C12E5', 'C1'): 'CH3', ('C12E5', 'C2'): 'CH2'}
    """
    dico_atoms = {}
    # Read the file that lists the atom
    with open(filename, 'r') as f :
        lines = f.readlines()
    # for each atom name (C1, C2 etc) of a given molecule (C12E5, BOG etc),
    # get the composition (CH3, CH2 etc)
    for l in lines:
        l = l.split()
        dico_atoms[(l[0],l[1])] = l[2]
    return dico_atoms


def asso_rad(dico_atom, dico_rad, table_atom):
    """
    Associate for each atom, its radius according to gromos

    -----------------------------------
    INPUT : 
    dico_atom : dictionnary
        Containing the group for each atom name of each molecule
    dico_rad : dictionnary
        Containing the gromos radius of the different groups
    table_atom : dataframe
        Containing every atom in the selection
    
    -----------------------------------
    OUTPUT
    dictionnary
        Containing for every atom its radius
        Example : {0: 0.2103, 1: 0.2284, etc}
    """
    dico_rad_num = {}
    for i in range(len(table_atom)):
        dico_rad_num[i] = dico_rad[dico_atom[(table_atom["resName"][i],table_atom["name"][i])]]
    return dico_rad_num


def find_index(coor, lim, size):
    """
    Find the index in the grid of a given position

    -----------------------------------
    INPUT : 
    coor : float
       The x or y coordinate of the studied atom minus or plus its radius
       to get its minimal or maximal position
    lim : float
        The x or y coordinate of the lowest atom
    size : float
        The column or row dimension of a cell
    
    -----------------------------------
    OUTPUT
    float
        The index in the grid of the extremity of an atom
    """
    # Compute the difference between the position coor and the minimal position
    # and then divide it by the size of a cell to find its index
    ind = np.floor((coor-lim)/size)
    return(ind)


def check_edge(ind, lim):
    """
    Check that the index is in the grid

    -----------------------------------
    INPUT : 
    ind : float
       The index in the grid (x or y) of an atom
    lim : float
        The index limit (in x or y)
    
    -----------------------------------
    OUTPUT
    float
        The index in the grid of the atom
    """
    # Check if the index is in the grid
    if lim == 0:
        # If out of the grid, return the first index in the grid
        if ind < lim:
            return lim
    else:
        # If out of the grid, return the last index in the grid
        if ind > lim:
            return lim
    return ind


def index_border(x, y, xmin, ymin, radius, colsize, rowsize):
    """
    Find the index of the outer atom (coordinate +- radius)
    and check that itt is in the grid

    -----------------------------------
    INPUT
    x : float
        The x coordinate of the studied atom
    y : float
        The y coordinate of the studied atom
    xmin : float
        The x coordinate of the lowest atom
    ymin : float
        The y coordinate of the lowest atom
    radius : float
        The radius of the studied atom
    colsize : float
        The column dimension of a cell
    rowsize : float
        The row dimension of a cell
    
    -----------------------------------
    OUTPUT
    tuple
        Containing the max and min indexes of an atom
    """
    # Find the index for the extremity of the atom radius
    col_min = find_index(x-radius, xmin, colsize)
    # Check that it is in the grid
    col_min = check_edge(col_min, 0.)
    col_max = find_index(x+radius, xmin, colsize)
    col_max = check_edge(col_max, args.column-1)
    row_min = find_index(y-radius, ymin, rowsize)
    row_min = check_edge(row_min, 0.)
    row_max = find_index(y+radius, ymin, rowsize)
    row_max = check_edge(row_max, args.row-1)
    return col_min, col_max, row_min, row_max


def search_pore_col(label_array, npore):
    """
    Search for the pore on the side of the matrix and label them as one
    if they are the same.

    This is due to the periodic conditions in MD. we observe the same pore on
    either side on the box so we need to label it as the same in the matrix

    -----------------------------------
    INPUT : 
    label_array : numpy matrix
       Containing the labeled pores
    npore : int
        Number of pores
    
    -----------------------------------
    OUTPUT
    numpy matrix
        Containing the reunited labeled pores
    int
        The new number of pores
    """
    # Get the list of row index where there is no atoms in the first column
    right =  np.where(label_array[:,0] > 0)[0]
    # Get the list of row index where there is no atoms in the last column
    left = np.where(label_array[:,-1] > 0)[0]
    
    # If there are more holes in the first column than the last one
    if len(right) >= len(left):
        # Loop on the index of the smallest list
        for index in left:
            # If the index also exists on the other side
            # and their label is different
            if index in right and label_array[index,-1] != label_array[index,0]:
                # Get the x,y index where this holes is labeled
                m = np.where(label_array==label_array[index,0])
                # Change the label on the right side for the one on the left
                # Basically, it is the same hole on both side
                # so we label it as one hole
                label_array[m] = label_array[index,-1]
                # Decrease the number of pores
                npore -= 1
    else:
        # Loop on the row index of the smallest list
        for index in right:
            # If the index also exists on the other side
            # and their label is different
            if index in left and label_array[index,-1] != label_array[index,0]:
                # Get the y,x index where this holes is labeled
                m = np.where(label_array==label_array[index,0])
                # Change the label on the right side for the one on the left
                label_array[m] = label_array[index,-1]
                # Decrease the number of pores
                npore -= 1
    return label_array, npore


def search_pore_row(label_array, npore):
    """
    Search for the pore on the side of the matrix and label them as one
    if they are the same.

    This is due to the periodic conditions in MD. we observe the same pore on
    either side on the box so we need to label it as the same in the matrix

    -----------------------------------
    INPUT : 
    label_array : numpy matrix
       Containing the labeled atom holes
    npore : int
        The number of pores
    
    -----------------------------------
    OUTPUT
    numpy matrix
        Containing the labeled pores
    int
        The new number of pores
    """
    # Get the list of column index where there is no atoms in the first row
    up =  np.where(label_array[0,] > 0)[0]
    # Get the list of column index where there is no atoms in the last row
    down = np.where(label_array[-1,] > 0)[0]

    # If there are more holes in the first row than the last one
    if len(up) >= len(down):
        # Loop on the index of the smallest list
        for index in down:
            # If the index also exists on the other side
            # and their label is different
            if index in up and label_array[-1, index] != label_array[0, index]:
                # Get the x,y index where this holes is labeled
                m = np.where(label_array==label_array[0,index])
                # Change the label on the right side for the one on the left
                # Basically, it is the same hole on both side
                # so we label it as one hole
                label_array[m] = label_array[-1,index]
                # Decrease the number of pores
                npore -= 1
    else:
        # Loop on the index of the smallest list
        for index in up:
            # If the index also exists on the other side
            # and their label is different
            if index in down and label_array[-1, index] != label_array[0, index]:
                # Get the x,y index where this holes is labeled
                m = np.where(label_array==label_array[0,index])
                # Change the label on the right side for the one on the left
                label_array[m] = label_array[-1,index]
                # Decrease the number of pores
                npore -= 1
    return label_array, npore


def get_labels(grid):
    """
    Label the space where there is no atoms and check that they are pores

    -----------------------------------
    INPUT
    grid : numpy matrix
       Containing the position of the pores
       0 : atoms
       1 : pores
    
    -----------------------------------
    OUTPUT
    numpy matrix
        Containing the reunited labeled pores
    int
        The number of pores
    """
    # Defining the structure of how the features connect
    # if none given, it will search for features that are in direct contact
    # Here, the features can also connect diagonally
    s = [[1,1,1],
         [1,1,1],
         [1,1,1]]
    # Labeling the grid (mlab) and counting the pores (npore)
    # Every non zero values in grid will be labeled
    mlab, npore = label(grid, structure=s)
    # Reunited the pores that are split bewteen the box edges
    mlab, npore = search_pore_col(mlab, npore)
    mlab, npore = search_pore_row(mlab, npore)

    return mlab, npore


def follow_pore(gridn, gridn1, labn, labn1):
    """
    Connect the actual pore labels to their previous one

    -----------------------------------
    INPUT
    gridn : numpy matrix
        Containing the previous grid matrix
        Contains 0 if first frame
    gridn1 : numpy matrix
        Containing the position of the pores
        0 : atoms
        1 : pores
    labn : numpy matrix
        Containing the previous label matrix
        Contains 0 if first frame
    labn1 : numpy matrix
        Containing the reuntied labeled pores
    
    -----------------------------------
    OUTPUT
    dictionnary
        Containing the actual labels as key and their previous one as value 
    """
    # Create a matrix adding the informations of gridn and gridn1
    grid = gridn + gridn1
    # get the rows and columns where the pore hasn't moved between two frames
    row,col = np.where(grid == 2)
    # Dictionnary that will allow us to know the different labels of
    # a same pore throughout the simulation
    connect = {}

    # Loop on the cells where the pore hasn't moved
    for i in range(len(row)):
        r = row[i]
        c = col[i]
        # If the label for the pore exists in connect
        if labn1[r,c] in connect:
            # If the previous value is higher than the previous value
            # This is in the case the actual pore regroups two previous ones
            if connect[labn1[r,c]] != labn[r,c] and connect[labn1[r,c]] > labn[r,c]:
                # Change to the lower previous value
                connect[labn1[r,c]] = labn[r,c]
        else:
            # The pore of label X now was the pore with the label Y previously
            connect[labn1[r,c]] = labn[r,c]
    return connect


def add_pore_in_time(connect, pore_in_time, frame):
    """
    Keep record of which frame was seen at which frame
    
    -----------------------------------
    INPUT
    connect : dictionnary
        Containing the actual labels as key and their previous one as value
    pore_in_time : dictionnary
        Containing the frames in which each pore labels is
    frame : int
        The number of the frame
    
    -----------------------------------
    OUTPUT
    dictionnary
        Containing the frames in which each pore labels is
        the labels are the keys and the lists of frames, the values
    """
    # Loop on the previous pore labels
    for pore in sorted(set(connect.values())):
        # If the previous label is in the dictionnary
        if pore in pore_in_time:
            # Add this frame to the pore label
            pore_in_time[pore] += [frame]
        else:
            # Add the label as key and the frame as value
            pore_in_time[pore] = [frame]
    return pore_in_time


def change_label(connect, maxlabel, labn1):
    """
    Change the actual pore labels to ttheir previous one

    This is  to keep the continuity of a pore
    
    -----------------------------------
    INPUT
    connect : dictionnary
        Containing the actual labels as key and their previous one as value
    maxlabel : int
        The highest label of the previous frame
    labn1 : numpy matrix
        Containing the reunited labeled pores
    
    -----------------------------------
    OUTPUT
    numpy matrix
        Containing the new labeled pore matrix
        with the pores having now the same label as their previous ones 
    Dictionnary
        Containing also the newly formed pores
    """
    # Get the rows and columns where there is a label 
    row,col = np.where(labn1 > 0)
    # Loop on these cells
    for i in range(len(row)):
        r = row[i]
        c = col[i]
        # If the label of a pore is in connect
        if labn1[r,c] in connect:
            # If the pore has divided in two
            if list(connect.values()).count(connect[labn1[r,c]]) > 1:
                # Don't change the label
                continue
            else:
                # Change the current label of a pore to its previous one
                labn1[r,c] = connect[labn1[r,c]]
        else:
            # Add this label to connect with a new label as value
            connect[labn1[r,c]] = maxlabel+1
            # Increase the maximum number of label because we just added one
            maxlabel += 1
            # Change the current label of a pore to its previous one
            labn1[r,c] = connect[labn1[r,c]]
    return labn1, connect


def add_dico(dico, key, value):
    """
    Add a value to a dictionnary with a specific key
    
    -----------------------------------
    INPUT
    dico : dictionnary
        Containing the surface of the pores
    key : int
        The key to be added
    value : float
        The value to be added
    
    -----------------------------------
    OUTPUT
    dictionnary
        Containing thge surface of the pores
        with the added value with the added key
    """
    if key in dico:
        dico[key] += value
    else:
        dico[key] = value
    return dico


def get_pores(labn, gridn, colsize, rowsize, indlabel, sizes):
    """
    Compute the surface of the pores and add it to a dictionnary
    
    -----------------------------------
    INPUT
    labn : numpy matrix
        Containing the reunited changed labeled pores of the previous frame
    gridn : numpy matrix
        Containing the position of the previous pores
        0 : atoms
        1 : pores
    colsize : float
        The column dimension of a cell of the previous frame
    rowsize : float
        The row dimension of a cell of the previous frame
    indlabel : list
        The list of the pore labels of the previous frame
    sizes : dictionnary
        Containing the real surface of the previous pores
    
    -----------------------------------
    OUTPUT
    list
        Containing the pore surfaces
    dictionnary
        Containing the real surface of the pores
    """
    # Compute the number of cells in grid that have been
    # labeled (in labels) as in indlabel 
    pores = sum(gridn, labn, indlabel)
    # Compute the surface of each pore
    pores = colsize*rowsize*pores
    # Loop on the pore labels
    for i, lab in enumerate(indlabel):
        # Add the surface to the sizes dictionnary
        sizes = add_dico(sizes, lab, [str(pores[i])])
    return pores, sizes


def write_pore_sizes(lpores, labpore, nframe):
    """
    Write a file containing the surface of each pore in a frame
    
    -----------------------------------
    INPUT
    lpores : list
        Containing the surface of the pores
    labpore : list
        Containing the pore labels
    nframe : int
        The number of the frame
    """
    with open("pores_size.txt", 'a') as f:
        f.write(f"FRAME {nframe}\n")
        f.write(f"number of pores:{len(lpores):5d}\n")
        for i, pore in enumerate(lpores):
            f.write(f"Pore number {labpore[i]}: {float(pore):.2f} nm^2\n")


def write_pdbs(nrow, ncol, csize, rsize, labels, xmin, ymin, zmax, nmdl):
    """
    Write pdb files containing the position of the pores
    
    -----------------------------------
    INPUT
    nrow : int
        The number of rows in the matrix
    ncol : int
        The number of columns in the matrix
    csize : float
        The column dimension of a cell
    rsize : float
        The column dimension of a cell
    labels : numpy matrix
        Containing the reunited changed labeled pores
    xmin : float
        The x coordinate of the lowest atom
    ymin : float
        The y coordinate of the lowest atom
    zmax : float
        The z of the highest atom at its highest during a simulation
    nmdl : int
        The number of the frame
    """
    # Open file to write inside
    with open(f"pdb/{str(nmdl)}.pdb","w") as f:
        f.write(f"MODEL {nmdl}\n")
        # Go through the label matrix 
        for i in range(int(nrow)):
            for j in range(int(ncol)):
                # If there is a pore
                if labels[i,j] > 0:
                    # convert the indexes of the pore to real number
                    x = xmin + j*csize + csize/2
                    y = ymin + i*rsize + rsize/2
                    # Write the poit as hetero atom in pdb format
                    f.write(f"HETATM{labels[i,j]:5d}  O   GRD A{labels[i,j]:4d}    {x*10:8.3f}{y*10:8.3f}{zmax*10:8.3f}\n")
        f.write("ENDMDL\n")


def apply_pbc(row, col):
    """
    Apply the periodic boundary condition to a pore in a frame

    -----------------------------------
    INPUT
    row : list
        Contains the rows in the label matrix of a pore
    col : int
        Contains the columns in the label matrix of a pore
    
    -----------------------------------
    OUTPUT
    list
        Containing the new rows in the label matrix of a pore
    list
	Containing the new columns in the label matrix of a pore
    """
    if max(row)-min(row) == args.row-1:
        if len(row[row < args.row/2]) < len(row[row > args.row/2]):
            # need to apply pbc
            row[np.where(row < args.row/2)] = row[row < args.row/2] + args.row
        else:
            row[np.where(row > args.row/2)] = row[row > args.row/2] - args.row        
    if max(col)-min(col) == args.column-1:
        if len(col[col < args.column/2]) < len(col[col > args.column/2]):
            # need to apply pbc
            col[np.where(col < args.column/2)] = col[col < args.column/2] + args.column
        else:
            col[np.where(col > args.column/2)] = col[col > args.column/2] - args.column
    
    return row, col


def get_com(labels, label_list, csize, rsize, xmin, ymin):
    """
    Compute the center of mass (COM) of each pore in a frame

    -----------------------------------
    INPUT
    labels : numpy matrix
        Containing the label matrix
        Contains 0 if first frame
    label_list : list
        Containing the pore labels
    csize : float
        The column dimension of a cell
    rsize : float
        The column dimension of a cell
    xmin : float
        The x coordinate of the lowest atom
    ymin : float
        The y coordinate of the lowest atom
    
    -----------------------------------
    OUTPUT
    dictionnary
        Containing the COM of the pores 
    """
    com = {}
    # Loop on the pore labels
    for label in label_list:
        # Get the rows and colummns where this pore is
        row,col = np.where(labels == label)
        # If the pore is across the simulation box
        # Apply PBC to the smallest part of the pore
        # And change the row and col lists
        row, col = apply_pbc(row, col)
        # Compute the mean column and row of this pore == COM of the pore
        # And convert the indexes into real values
        ymean = xmin + np.mean(row) * rsize + rsize/2
        xmean = ymin + np.mean(col) * csize + csize/2
        # Add the COM of the pore to a dictionary
        com[label] = [xmean, ymean]
    return com


def write_com(com, zmax, frame, xmax, ymax):
    """
    Write the center of mass (COM) of each pore in a frame in a file

    -----------------------------------
    INPUT
    com : dictionnary
        Containing the COM of the pores in a frame
    zmax : float
        The z of the highest atom at its highest during a simulation
    frame : int
        The number of the frame
    """
    # Open file to write inside
    with open("center_of_mass_pore.pdb", 'a') as f:
        # Write first line in pdb to apply pbc in vmd
        f.write(f"TITLE     Pore system t=   {(frame-1)*1000}.00000 step= {args.step*(frame-1)}\n")
        f.write(f"CRYST1  {xmax*10:.3f}  {ymax*10:.3f}   {zmax*10:.3f}  90.00  90.00  90.00 P 1           1\n")
        f.write(f"MODEL {frame}\n")
        # Get the pore labels sorted
        keys = sorted(com.keys())
        # Write in the file for each pore label
        for k in keys:
            f.write(f"HETATM{k:5d}  O   GRD A{k:4d}    {com[k][0]*10:8.3f}{com[k][1]*10:8.3f}{zmax*10:8.3f}\n")
        f.write("ENDMDL\n")


def write_empty_pdbs(nmdl):
    """
    Write a pdb file containing nothing
    
    -----------------------------------
    INPUT
    nmdl : int
        The number of the frame
    """ 
    # Open file to write inside
    with open(f"pdb/{str(nmdl)}.pdb","a") as f:
        f.write(f"MODEL {nmdl}\n")
        f.write("ENDMDL\n")


def write_stat_pore(dico_pore, steps):
    """
    Write a file containing the first and last frame of each pore
    Each line is : 
        pore label;first frame;last frame;number of frame;percentage of frame with this pore
    
    -----------------------------------
    INPUT
    dico_pore : dictionnary
        Containing the frames in which each pore labels is
    steps : int
        The number of frames
    """
    # Open file to write inside
    with open("pores_in_time.dat", 'w') as f:
        # Get all the sorted pore labels
        keys = sorted(dico_pore.keys())
        # Write in the file for each pore label
        for k in keys:
            f.write(f"{k};{dico_pore[k][0]};{dico_pore[k][-1]};{len(dico_pore[k])};{len(dico_pore[k])/float(steps)*100:.2f}\n")


def write_mat_size(dico_pore, sizes, frames):
    """
    Write a file containing the matrix of pore size in nm^2 for each frame.
        Lines: pore (first line is bilayer area in nm^2)
        Columns: frame (first column is pore name)
    
    -----------------------------------
    INPUT
    dico_pore : dictionnary
        Containing the frames in which each pore labels is
    sizes : dictionnary
        Containing the real surface of the pores
    frames : int
        The number of frames
    """
    # Open file to write inside
    with open("sizes.dat", 'w') as f:
        # Write the bilayer area
        f.write(f"0;{sizes['0'][0]}\n")
        # Get the sorted pore labels
        pores = sorted(dico_pore.keys())
        # Loop on the labels
        for p in pores:
            # If there is a surface for this pore
            if p in sizes:
                # Write its area in the file
                line = ["0"]*(frames+1)
                line[0] = str(p)
                line[dico_pore[p][0]:dico_pore[p][-1]+1] = sizes[p]
                line = ";".join(line)
                f.write(line+"\n")


################################################################################
##################################### MAIN #####################################
################################################################################

if __name__=="__main__":
    gromos_radius = {"CH3":0.2103, "CH2":0.2284, "OE":0.1871, "OA":0.1744, "CH1r":0.2817, "Or":0.1599}
    directory = "pdb"

    # Get arguments
    args = recup_args()

    # Create directory for pdb files
    make_dir(directory)

    # Load trajectory
    print("Trajectory is loading...")
    if args.topology:
        traj = md.load(args.trajectory, top=args.topology)
    else:
        traj = md.load(args.trajectory)

    # Select in atoms in the index file
    selec = read_ndx(args.index)
    traj = traj.atom_slice(selec)

    # Get the topology in the form of a dataframe
    topo = traj.topology
    topo_dtf = topo.to_dataframe()[0]

    # Get the molecule and general name of each atom group
    dico_atom = get_atoms(args.list_atoms)
    # Get the radius of each atom group
    dico_rad = asso_rad(dico_atom, gromos_radius, topo_dtf)

    # Get all the coordinates (x,y,z) of all the atoms for each frame
    coor = traj.xyz
    # Get the z coordinates of the highest atom at its highest in the simulation
    z = coor.max(1).max(0)[2]

    # Initialise dictionnaries
    pore_in_time = {}
    sizes = {"0":[]}
    # Initialise two matrixes row*column at 0
    # one for the labeled pores of the previous frame
    # and one for the grid containing the pores of the previous frame
    lab_prev = np.zeros((int(args.row),int(args.column)), dtype=int)
    grid_prev = np.zeros((int(args.row),int(args.column)), dtype=int)

    # Get the x and y coordinates of the highest and lowest atom
    xmin_prev, ymin_prev = coor[0].min(axis=0)[:2]  # 0 colum 1 row
    xmax_prev, ymax_prev = coor[0].max(axis=0)[:2]
    # Defining the cell size    
    csize_prev = (xmax_prev-xmin_prev)/args.column
    rsize_prev = (ymax_prev-ymin_prev)/args.row
    
    # Write first line in pdb to apply pbc in vmd
    with open("center_of_mass_pore.pdb", 'w') as f :
        f.write("CRYST1  136.037  136.037   94.844  90.00  90.00  90.00 P 1           1\n")

    # Loop on the frames
    for frame, crd in enumerate(coor):
        sys.stdout.write(f" FRAME {frame+1:5d}/{len(coor):5d}\r")
        sys.stdout.flush()

        # Get the x and y coordinates of the highest and lowest atom
        xmin, ymin = crd.min(axis=0)[:2]  # 0 colum 1 row
        xmax, ymax = crd.max(axis=0)[:2]
        # Initialise a matrix row*column at 1
        grid = np.ones((int(args.row),int(args.column)), dtype=int)
        # Defining the cell size    
        csize = (xmax-xmin)/args.column
        rsize = (ymax-ymin)/args.row
        # Add the bilayer area in real number at this frame
        # to the sizes dictionnary
        sizes["0"] += [str(csize*rsize*args.column*args.row)]
        
        # Loop on the atoms
        for i,atm in enumerate(crd):
            # Get the radius of this atom
            rad = dico_rad[i]
            # Get the indexes of the extremities of the atom
            c1, c2, r1, r2 = index_border(atm[0], atm[1], xmin, ymin, rad, csize, rsize)
            # Fill the matrix with 0 where the atom is
            grid[int(r1):int(r2+1), int(c1):int(c2+1)] = 0
        
        # Label the pores, conect them across the simulation box ans count them
        labels, npore = get_labels(grid)

        # Connect the actual pore labels to their previous ones : eliminates
        # the small ones
        connect = follow_pore(grid_prev, grid, lab_prev, labels)
        # Log which pore labels were seen at which frame
        pore_in_time = add_pore_in_time(connect, pore_in_time, frame)

        # if pore_in_time isn't empty
        if pore_in_time:
            # Change the actual pore labels to their previous ones expect if
            # one pore has dived in two
            labels, connect = change_label(connect, max(pore_in_time), labels)
         
        # Get the list of the pore labels for the previous frame
        list_of_pore = list(set(pore_in_time.keys())&set(connect.values()))
        # Add the surface of the pores to sizes dictionnary for the previous frame
        pores, sizes = get_pores(lab_prev, grid_prev, csize_prev, rsize_prev, list_of_pore, sizes)
        # Write a .txt file to keep trace of the pore surfaces
        write_pore_sizes(pores, list_of_pore, frame)
        
        # Write pdbs with the position of the previous pores
        if pores.any():
            write_pdbs(args.row, args.column, csize_prev, rsize_prev, lab_prev, xmin_prev, ymin_prev, z, frame)
            # Write the COM of the pores in a .txt file
            com = get_com(lab_prev, list_of_pore, csize_prev, rsize_prev, xmin_prev, ymin_prev)
            write_com(com, z, frame, xmax_prev, ymax_prev)
        else:
            write_empty_pdbs(frame)

        # Update the previous label, grid matrix etc
        lab_prev = labels
        grid_prev = grid
        csize_prev = csize
        rsize_prev = rsize
        xmin_prev = xmin
        ymin_prev = ymin
        xmax_prev = xmax
        ymax_prev = ymax
    
    write_stat_pore(pore_in_time, len(coor))
    write_mat_size(pore_in_time, sizes, len(coor))
    print("")
