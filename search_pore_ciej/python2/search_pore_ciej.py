#!/usr/bin/python
# -*- coding: utf-8 -*-

__author__ = "Caroline SENAC"
__date__ = "10/2016"
__version__ = "V3"

import re
import argparse
import sys
import os
import errno

from scipy.ndimage import measurements
import numpy as np
import mdtraj as md
import pandas as pd

usageprog = \
"""
    python search_pore_ciej.py -f traj.xtc -c topol.gro -n index.ndx -p list_atoms.txt [-col 200] [-row 200]
"""


def recup_args():
    parser = argparse.ArgumentParser(description="Find pore in CiEj bilayer", usage=usageprog)
    parser.add_argument("-f", "--trajectory", help="<.xtc/.gro> trajectory file")
    parser.add_argument("-c", "--topology", help="<.gro> topology file")
    parser.add_argument("-n", "--index", help="<.ndx> index file")
    parser.add_argument("-p", "--list_atoms", help="<.txt> list of atoms")
    parser.add_argument("-col", "--column", type=float, default=200.,\
                        help="<int> number of column in grid")
    parser.add_argument("-row", "--row", type=float, default=200.,\
                        help="<int> number of row in grid")
    args = parser.parse_args()
    return args


def make_dir(dirpath):
    try:
        os.makedirs(dirpath)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
        else:
            print "Directory {:s} already exists!".format(dirpath)


def read_ndx(filename):
    l_index = []
    regex = re.compile("\[")
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()
    for l in lines:
        if not regex.match(l):
            l_index += [int(i) for i in l.split()]
    l_index = np.array(sorted(l_index)) -1
    return l_index


def recup_atoms(filename):
    dico_atoms = {}
    f = open(filename, 'r')
    lines = f.readlines()
    for l in lines:
        l = l.split()
        dico_atoms[(l[0],l[1])] = l[2]
    return dico_atoms


def asso_rad(dico_atom, dico_rad, table_atom):
    dico_rad_num = {}
    for i in xrange(len(table_atom)):
        dico_rad_num[i] = dico_rad[dico_atom[(table_atom["resName"][i],table_atom["name"][i])]]
    return dico_rad_num


def find_index(coor, lim, size):
    ind = np.floor((coor-lim)/size)
    return(ind)


def check_edge(ind, lim):
    if lim == 0:
        if ind < lim:
            return lim
    else:
        if ind > lim:
            return lim
    return ind


def search_pore_row(label_array):
    up =  np.where(label_array[0,] > 0)[0]
    down = np.where(label_array[-1,] > 0)[0]
    if len(up) >= len(down):
        for e in down:
            if e in up:
                m = np.where(label_array==label_array[0,e])
                label_array[m] = label_array[-1,e]
    else:
        for e in up:
            if e in down:
                m = np.where(label_array==label_array[0,e])
                label_array[m] = label_array[-1,e]
    return label_array


def search_pore_col(label_array):
    up =  np.where(label_array[:,0] > 0)[0]
    down = np.where(label_array[:,-1] > 0)[0]
    if len(up) >= len(down):
        for e in down:
            if e in up:
                m = np.where(label_array==label_array[e,0])
                label_array[m] = label_array[e,-1]
    else:
        for e in up:
            if e in down:
                m = np.where(label_array==label_array[e,0])
                label_array[m] = label_array[e,-1]
    return label_array


def get_labels(grid):
    s = [[1,1,1],
         [1,1,1],
         [1,1,1]]
    mlab, npore = measurements.label(grid, structure=s)
    mlab = search_pore_col(mlab)
    mlab = search_pore_row(mlab)
    return mlab, npore


def index_border(x, y, xmin, ymin, radius, colsize, rowsize):
    col_min = find_index(x-radius, xmin, colsize)
    col_min = check_edge(col_min, 0.)
    col_max = find_index(x+radius, xmin, colsize)
    col_max = check_edge(col_max, args.column-1)
    row_min = find_index(y-radius, ymin, rowsize)
    row_min = check_edge(row_min, 0.)
    row_max = find_index(y+radius, ymin, rowsize)
    row_max = check_edge(row_max, args.row-1)
    return col_min, col_max, row_min, row_max


def follow_pore(gridn, gridn1, labn, labn1):
    grid = gridn + gridn1
    row,col = np.where(grid == 2)
    connect = {}
    for i in xrange(len(row)):
        r = row[i]
        c = col[i]
        if labn1[r,c] in connect:
            if connect[labn1[r,c]] > labn[r,c]:
                connect[labn1[r,c]] = labn[r,c]
        else:
            connect[labn1[r,c]] = labn[r,c]
    return connect


def add_pore_in_time(connect,pore_in_time, frame):
    for pore in sorted(set(connect.values())):
        if pore in pore_in_time:
            pore_in_time[pore] += [frame]
        else:
            pore_in_time[pore] = [frame]
    return pore_in_time


def change_label(connect, maxlabel, labn1):
    row,col = np.where(labn1 > 0)
    for i in xrange(len(row)):
        r = row[i]
        c = col[i]
        if labn1[r,c] in connect:
            labn1[r,c] = connect[labn1[r,c]]
        else:
            connect[labn1[r,c]] = maxlabel+1
            maxlabel += 1
            labn1[r,c] = connect[labn1[r,c]]
    return labn1, connect


def get_pores(labels, grid, colsize, rowsize, npore, indlabel, sizes):
    pores = measurements.sum(grid, labels, indlabel)
    pores = pores*colsize*rowsize
    for i,k in enumerate(indlabel):
        sizes = add_dico(sizes, k, [str(pores[i])])
    return pores, sizes


def add_dico(dico, k, v):
    if k in dico:
        dico[k] += v
    else:
        dico[k] = v
    return dico


def print_mat_size(times, sizes, frames):
    f = open("sizes.dat", 'w')
    line = "0;"
    line = line+(";".join(sizes["0"]))
    f.write(line+"\n")
    pores = sorted(times.keys())
    for p in pores:
        if p in sizes:
            line = ["0"]*(frames+1)
            line[0] = str(p)
            line[times[p][0]:times[p][-1]+1] = sizes[p]
            line = ";".join(line)
            f.write(line+"\n")
    f.close()


def print_output(lpores, labpore, nframe):
    f = open("pores_size.txt", 'a')
    s = "FRAME {:d}\n".format(nframe)
    f.write(s)
    f.write("number of pores:{:5d}\n".format(len(lpores)))
    for p, pore in enumerate(lpores):
        f.write("Pore number {:d}: {:.2f} nm^2\n".format(labpore[p], float(pore)))
    f.close()


def write_pdb(nrow, ncol, csize, rsize, labels, xmin, ymin, zmax, nmdl):
    f = open("pores.pdb","a")
    s = "MODEL {:d}\n".format(nmdl)
    f.write(s)
    for i in xrange(int(nrow)):
        for j in xrange(int(ncol)):
            if labels[i,j] > 0:
                x = xmin + j*csize + csize/2
                y = ymin + i*rsize + rsize/2
                s = "HETATM{:5d}  O   GRD A{:4d}    {:8.3f}{:8.3f}{:8.3f}\n".format(labels[i,j],labels[i,j],x*10,y*10,zmax*10)
                f.write(s)
    f.write("ENDMDL\n")
    f.close()


def write_pdbs(nrow, ncol, csize, rsize, labels, xmin, ymin, zmax, nmdl):
    name="pdb/"+str(nmdl)+".pdb" 
    f = open(name,"w")
    s = "MODEL {:d}\n".format(nmdl)
    f.write(s)
    for i in xrange(int(nrow)):
        for j in xrange(int(ncol)):
            if labels[i,j] > 0:
                x = xmin + j*csize + csize/2
                y = ymin + i*rsize + rsize/2
                s = "HETATM{:5d}  O   GRD A{:4d}    {:8.3f}{:8.3f}{:8.3f}\n".format(labels[i,j],labels[i,j],x*10,y*10,zmax*10)
                f.write(s)
    f.write("ENDMDL\n")
    f.close()


def write_empty_pdbs(nmdl):
    name="pdb/"+str(nmdl)+".pdb"   
    f = open(name,"a")
    s = "MODEL {:d}\n".format(nmdl)
    f.write(s)
    f.write("ENDMDL\n")
    f.close()


def write_trou(dico_trou, steps):
    f = open("pores_in_time.dat", 'w')
    key = sorted(dico_trou.keys())
    for k in key:
        f.write("{:d};{:d};{:d};{:d};{:.2f}\n".format(k,dico_trou[k][0],\
                                           dico_trou[k][-1],len(dico_trou[k]),\
                                           len(dico_trou[k])/float(steps)*100))
    f.close()


################################################################################
##################################### MAIN #####################################
################################################################################

gromos_radius = {"CH3":0.2103, "CH2":0.2284, "OE":0.1871, "OA":0.1744, "CH1r":0.2817, "Or":0.1599}
directory = "pdb"

#recup args
args = recup_args()

#create directory for pdb files
make_dir(directory)

#load trajectory
print "Trajectory is loading..."
if args.topology:
    traj = md.load(args.trajectory, top=args.topology)
else:
    traj = md.load(args.trajectory)

#selection
selec = read_ndx(args.index)
traj = traj.atom_slice(selec)

topo = traj.topology
table = topo.to_dataframe()[0]

dico_at = recup_atoms(args.list_atoms)
dico_rad = asso_rad(dico_at, gromos_radius, table)

coor = traj.xyz
z = coor.max(1).max(0)[2]
pore_in_time = {}
sizes = {"0":[]}
labn = np.zeros((int(args.row),int(args.column)), dtype=int)
gridn = np.zeros((int(args.row),int(args.column)), dtype=int)

for frame, crd in enumerate(coor):
    status = " FRAME {:5d}/{:5d}".format(frame+1, len(coor))
    #print frame+1
    sys.stdout.write('{:s}\r'.format(status))
    sys.stdout.flush()
    #searching min and max surface
    xmin, ymin = crd.min(axis=0)[:2]  # 0 colum 1 row
    xmax, ymax = crd.max(axis=0)[:2]
    grid = np.ones((int(args.row),int(args.column)), dtype=int)
    #cell size    
    csize = (xmax-xmin)/args.column
    rsize = (ymax-ymin)/args.row
    #sizes    
    sizes["0"] += [str(csize*rsize*args.column*args.row)]
    for i,atm in enumerate(crd):
        rad = dico_rad[i]
        c1, c2, r1, r2 = index_border(atm[0], atm[1], xmin, ymin, rad, csize, rsize)
        grid[int(r1):int(r2+1), int(c1):int(c2+1)] = 0
    
    labels, npore = get_labels(grid)
    connect = follow_pore(gridn, grid, labn, labels)
    pore_in_time = add_pore_in_time(connect,pore_in_time, frame+1)
    if pore_in_time:
        labels, connect = change_label(connect, max(pore_in_time), labels)
    list_of_pore = list(set(pore_in_time.keys())&set(connect.values()))
    pores, sizes = get_pores(labels, grid, csize, rsize, npore, list_of_pore, sizes)
    print_output(pores, list_of_pore, frame+1)
    if pores.any():
        write_pdbs(args.row, args.column, csize, rsize, labels, xmin, ymin, z, frame+1)
    else:
        write_empty_pdbs(frame+1)
    labn = labels
    gridn = grid

write_trou(pore_in_time, len(coor))
print_mat_size(pore_in_time, sizes, len(coor))
print ""


