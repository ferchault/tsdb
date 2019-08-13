import sys
import numpy as np


def get_mol(filename):
  ''' read in labels and coordinates
  '''
  lines = open(filename, 'r').readlines()

  labels = []
  coords = []

  for i, line in enumerate(lines):
    if i == 0: numAtoms = line
    if i == 1: continue

    if i > 1:
      tokens = line.split()
      labels.append(tokens[0])
      coords.append([ float(tokens[1]), float(tokens[2]), float(tokens[3]) ])

  return numAtoms, labels, coords


def print_xyz(numAtoms, labels, coords):
  ''' print xyz file
  '''
  print(numAtoms)

  for i in range(len(coords)):
    print(labels[i], coords[i][0], coords[i][1], coords[i][2])

def is_ok_CX_vs_XH(labels, coords):
  corrupted_XH_mols = []
  C1 = coords[0]
  X  = coords[2]
  CX = np.array([ C1[0] - X[0], C1[1] - X[1], C1[2] - X[2] ])
  CX_len = np.linalg.norm(CX)

  for i in range(len(labels)):
    if labels[i] == "H":
      H = coords[i]
      XH = np.array([ H[0] - X[0], H[1] - X[1], H[2] - X[2] ])
      XH_len = np.linalg.norm(XH)

      if CX_len > XH_len: return False

  return True

if __name__ == "__main__":

  filename =  sys.argv[1]

  numAtoms, labels, coords = get_mol(filename)

  test = is_ok_CX_vs_XH(labels, coords)

  sys.exit(1-int(test))





