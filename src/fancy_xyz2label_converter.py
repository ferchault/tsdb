import sys
import math
#import pybel
import numpy as np
import openbabel as ob

def translate(xyz):
	xyz -= xyz[0]

	return xyz

def rotate(xyz, axis, theta):
	rotM = rotation_matrix(axis, theta)

	for i in range(1,len(xyz)):
		xyz[i] = np.dot(rotM, xyz[i]) 

	return xyz

def rotation_matrix(axis, theta):
	""" Return the rotation matrix associated with counterclockwise rotation about
			the given axis by theta radians.
	"""
	axis = axis/np.linalg.norm(axis)

	a = np.cos(theta/2.0)
	b, c, d = -axis*np.sin(theta/2.0)
	aa, bb, cc, dd = a*a, b*b, c*c, d*d
	bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
	
	return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                 [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                 [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

def print_geom(xyz):
	print("\n")
	for i,atom in enumerate(xyz):
 		print label[i], atom[0], atom[1],atom[2]

def write_geom(xyz):
	f = open("tmp.xyz", 'w')
	f.write(str(nAtoms) + "\n")
	f.write('\n')
	for i, atom in enumerate(xyz):
		string = str(label[i] + "\t" + str(atom[0]) + "\t" +  str(atom[1]) + "\t" + str(atom[2]) + "\n")
		f.write(string)
	f.close()

def get_norm(vec):
	return vec / np.linalg.norm(vec)

def get_theta(vec1, vec2):
	return -np.arccos(np.dot(vec1, vec2))

def get_axis(vec1, vec2):
	return np.cross(vec1, vec2)

def check_FG(mol, nbrAtom):
	Atoms = ob.OBAtom
	count = 0
	for neighbor in ob.OBAtomAtomIter(nbrAtom):
		count += 1

	if count == 1:
		return "A"
	if count == 2:
		return "C"
	if count == 4:
		return "D"

	if count == 3:
		fg = "B"
		for neighbor in ob.OBAtomAtomIter(nbrAtom):
			if mol.GetAtom(neighbor.GetIdx()).GetAtomicNum() == 1:
				fg = "E"
		return fg 

if __name__ == "__main__":

	# read in data
	filename = sys.argv[1]

	f = open(filename, 'r')
	lines = f.readlines()
	f.close()

	# declare arrays
	label = np.asarray([])
	xyz = np.asarray([])
	ex = np.asarray([1,0,0])
	ez = np.asarray([0,0,1])

	for i,line in enumerate(lines):
		if i == 0: nAtoms = int(line)
		if i == 1: continue
		if i == 4: X = line.split()[0]
		if i == 5: Y = line.split()[0]
		if i > 1:
			tokens = line.split()
			label  = np.append(label, tokens[0])
			xyz    = np.append(xyz, [float(tokens[1]), float(tokens[2]), float(tokens[3])])

	xyz = np.reshape(xyz,(nAtoms,3))

	######################
	# translate and rotate
	######################

	xyz = translate(xyz)

	CC = xyz[1]-xyz[0]
	eCC = get_norm(CC)
	theta_C = get_theta(ex, eCC)
	axis_C  = get_axis(ex, CC)

	xyz = rotate(xyz, axis_C, theta_C)

	CX = xyz[2]-xyz[0]
	Vec = [0, CX[1], CX[2]]
	eVec = get_norm(Vec)
	theta_X = get_theta(ez, eVec)
	axis_X  = ex

	xyz = rotate(xyz, axis_X, theta_X)

	#print_geom(xyz)
	write_geom(xyz)

	######################
	# open babel stuff 
	######################

	mol = ob.OBMol()
	conv = ob.OBConversion()
	conv.SetInAndOutFormats("xyz","xyz")
	conv.ReadFile(mol, "tmp.xyz")
	C1 = ob.OBAtom
	C2 = ob.OBAtom
	nbrAtom = ob.OBAtom

	C1 = mol.GetAtom(1)
	C2 = mol.GetAtom(2)

	R1,R2,R3,R4 = ["0" for i in range(4)]
	for neighbor in ob.OBAtomAtomIter(C1):
		nbr = neighbor.GetIdx()
		nbrAtom = mol.GetAtom(nbr)
		if nbr == 2: continue
		if nbr == 3: continue
		vec = [nbrAtom.GetX(), nbrAtom.GetY(), nbrAtom.GetZ()]
		if vec[1] > 0:
			R1 = check_FG(mol, nbrAtom)
		if vec[1] < 0:
			R2 = check_FG(mol, nbrAtom)

	for neighbor in ob.OBAtomAtomIter(C2):
		nbr = neighbor.GetIdx()
		nbrAtom = mol.GetAtom(nbr)
		if nbr == 1: continue
		if nbr == 5: continue
		vec = [nbrAtom.GetX(), nbrAtom.GetY(), nbrAtom.GetZ()]
		if vec[1] > 0:
			R4 = check_FG(mol, nbrAtom)
		if vec[1] < 0:
			R3 = check_FG(mol, nbrAtom)

	if X == "F": X = "A"
	if X == "Cl": X = "B"
	if X == "Br": X = "C"
	if Y == "H": Y = "A"
	if Y == "F": Y = "B"
	if Y == "Cl": Y = "C"
	if Y == "Br": Y = "D"

	print R1 + "_" + R2 + "_" + R3 + "_" + R4 + "_" + X + "_" + Y

