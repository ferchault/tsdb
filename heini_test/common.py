import sys
import time
from datetime import datetime
import random
import cPickle
import numpy as np
from copy import deepcopy
import qml
from qml.representations import *
from qml.kernels import gaussian_kernel
from qml.kernels import laplacian_kernel
from qml.fchl import get_local_symmetric_kernels
from qml.math import cho_solve
import itertools
from time import time
from qml.distance import l2_distance

'''
def get_FCHL(mols):
  for mol in mols:
    mol.representation = generate_representation(mol.coordinates,mol.nuclear_charges,max_size=123)

  X = np.asarray([mol.representation for mol in mols])
  Y = np.asarray([ mol.properties for mol in mols ])

  return X, Y
'''
def get_SLATM(mols):
  mbtypes = get_slatm_mbtypes([mol.nuclear_charges for mol in mols])
  for mol in mols:
    mol.generate_slatm(mbtypes)

  X = np.asarray([mol.representation for mol in mols])
  Y = np.asarray([ mol.properties for mol in mols ])

  return X, Y

def get_CM(mols):
  for mol in mols:
    mol.generate_coulomb_matrix(size=132)

  X = np.asarray([mol.representation for mol in mols])
  Y = np.asarray([ mol.properties for mol in mols ])

  return X, Y


def get_BoB(mols):
  start = time()

  bags = {
      "H":  max([mol.atomtypes.count("H" ) for mol in mols]),
      "C":  max([mol.atomtypes.count("C" ) for mol in mols]),
      "N":  max([mol.atomtypes.count("N" ) for mol in mols]),
      "O":  max([mol.atomtypes.count("O" ) for mol in mols]),
      "F":  max([mol.atomtypes.count("F" ) for mol in mols]),
      "Si": max([mol.atomtypes.count("Si") for mol in mols]),
      "Br": max([mol.atomtypes.count("Br") for mol in mols]),
      "P":  max([mol.atomtypes.count("P" ) for mol in mols]),
      "Cl": max([mol.atomtypes.count("Cl") for mol in mols]),
      "Ni": max([mol.atomtypes.count("Ni") for mol in mols]),
      "Pt": max([mol.atomtypes.count("Pt") for mol in mols]),
      "Pd": max([mol.atomtypes.count("Pd") for mol in mols]),
      "Cu": max([mol.atomtypes.count("Cu") for mol in mols]),
      "Ag": max([mol.atomtypes.count("Ag") for mol in mols]),
      "Au": max([mol.atomtypes.count("Au") for mol in mols])
    }

  for mol in mols:
    mol.generate_bob(asize=bags)

  end = time()
  print(end-start)

  X = np.asarray([mol.representation for mol in mols])
  Y = np.asarray([ mol.properties for mol in mols ])

  return X, Y

def get_properties(filename):
  print "\n -> get compound class"
  f = open(filename, "r")
  lines = f.readlines()
  f.close()

  properties = dict()

  for line in lines:
    tokens = line.split()
    name = tokens[0]
    property = float(tokens[1])
    properties[name] = property

  return properties

def get_molecules(mols_dict, data):
  mols = []
  for xyz_file in mols_dict:
    mol = qml.Compound()
    mol.read_xyz("xyz/" + xyz_file + ".xyz")
    mol.properties = data[xyz_file]
    mols.append(mol)
  total = len(mols)

  return mols, total

def get_Representation(mols, Rep):
  print "\n -> calculating the Representation"
  if Rep == "BoB": return get_BoB(mols)
  if Rep == "CM": return get_CM(mols)
  if Rep == "SLATM": return get_SLATM(mols)
  #if Rep == "FCHL": return get_FCHL(mols)
  else:
    print("\ninvalid representation\n")
    exit(1)

def get_Kernel(X, sigmas, kernel):
  print "\n -> calculating the Kernel"
  start = time()
  if kernel == "Gaussian": K = gaussian_kernel(X,X,sigmas)
  if kernel == "Laplacian": K = laplacian_kernel(X,X,sigmas)
  #if kernel == "local_symmetric" : K = get_local_symmetric_kernels(X, kernel_args={"sigma":sigmas}, alchemy="off", two_body_power=4.0, three_body_power=2.0, cut_start=0.6, cut_distance=6.0)
  else:
    print("\ninvalid Kernel\n")
    exit(1)
  end = time()
  print(end-start)

  return K
'''
def CrossValidation_fchl(K,train, nModels, llambda,j, total, Yprime):
  total = len(parameters["total"])
  test = total - train
  maes = []

  for i in range(nModels):
    split = range(total)
    random.shuffle(split)

    training_index  = split[:train]
    test_index      = split[-test:]

    Y = Yprime[training_index]
    Ys = Yprime[test_index]

    C = deepcopy(K[j][training_index][:,training_index])
    C[np.diag_indices_from(C)] += llambda

    alpha = cho_solve(C, Y)

    Yss = np.dot((K[j][training_index][:,test_index]).T, alpha)
    diff = Yss  - Ys
    mae = np.mean(np.abs(diff))
    maes.append(mae)
  s = np.std(maes)/np.sqrt(nModels)

  return maes, s
'''
def CrossValidation(K,train, nModels, llambda, total, Yprime):
  test = total - train
  maes = []

  for i in range(nModels):
    split = range(total)
    random.shuffle(split)

    training_index  = split[:train]
    test_index      = split[-test:]

    Y = Yprime[training_index]
    Ys = Yprime[test_index]

    C = deepcopy(K[training_index][:,training_index])
    C[np.diag_indices_from(C)] += llambda

    alpha = cho_solve(C, Y)

    Yss = np.dot((K[training_index][:,test_index]).T, alpha)
    diff = Yss  - Ys
    mae = np.mean(np.abs(diff))
    maes.append(mae)
  s = np.std(maes)/np.sqrt(nModels)

  return maes, s
