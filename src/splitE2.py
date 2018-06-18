#!/usr/bin/env python

import sys
import os
import numpy as np
#import cclib

from cclib.io import ccopen

#from cclib.parser import ccopen

NAME = {
    1:  'H'  ,
    2:  'He' ,
    3:  'Li' ,
    4:  'Be' ,
    5:  'B'  ,
    6:  'C'  ,
    7:  'N'  ,
    8:  'O'  ,
    9:  'F'  ,
   10:  'Ne' ,
   11:  'Na' ,
   12:  'Mg' ,
   13:  'Al' ,
   14:  'Si' ,
   15:  'P'  ,
   16:  'S'  ,
   17:  'Cl' ,
   18:  'Ar' ,
   19:  'K'  ,
   20:  'Ca' ,
   21:  'Sc' ,
   22:  'Ti' ,
   23:  'V'  ,
   24:  'Cr' ,
   25:  'Mn' ,
   26:  'Fe' ,
   27:  'Co' ,
   28:  'Ni' ,
   29:  'Cu' ,
   30:  'Zn' ,
   31:  'Ga' ,
   32:  'Ge' ,
   33:  'As' ,
   34:  'Se' ,
   35:  'Br' }


if __name__ == "__main__":
    # read in log files
    filename = sys.argv[1]
    #filename = "run.log"
    mylogfile = ccopen(filename)
    data = mylogfile.parse()
    
    # extract data
    coords = data.atomcoords
    nAtoms = data.natom
    atoms = [NAME[i] for i in data.atomnos]
    
    # take last geometry
    coords = coords[-1]

    
    if np.linalg.norm(coords[1] - coords[4]) < 1.5:
        coords = np.delete(coords, 3, 0)
        atoms = np.delete(atoms, 3, 0)
        side = "reactant"
    else:
        coords = np.delete(coords, 4, 0)
        coords = np.delete(coords, 3, 0)
        coords = np.delete(coords, 2, 0)
        atoms = np.delete(atoms, 4, 0)
        atoms = np.delete(atoms, 3, 0)
        atoms = np.delete(atoms, 2, 0)
        side = "product"     
    
    if side == "product":
        print(nAtoms - 3)
        print(side, "1")
        for i,atom in enumerate(coords):
            print(str(atoms[i]), *atom)
    
    if side == "reactant":
        print(nAtoms - 1)
        print(side, "-1")
        for i,atom in enumerate(coords):
            print(str(atoms[i]), *atom)
