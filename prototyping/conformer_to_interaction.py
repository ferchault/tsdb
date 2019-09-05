#!/bin/bash
# Usage: conformer_to_interaction.py XYZFILE
# e.g.   conformer_to_interaction.py run.xyz
# builds the following input structure in the current directory for interaction energy optimisations:
# int-RXN-YLABEL

import sys
import os
import numpy as np
from scipy.spatial.transform import Rotation

def read_xyz(fn):
	with open(fn) as fh:
		lines = fh.readlines()[2:]
		elements = []
		coords = []
	for line in lines:
		parts = line.strip().split()
		elements.append(parts[0])
		coords.append([float(_) for _ in parts[1:])

	return elements, np.array(coords)

def do_it(elem, coord, rxn, ylabel):
	if rxn == 'e2':
		base_idx = 1
		pos_idx = 3
		ref_idx = 0
		addX = 2.3
	if rxn == 'sn2':
		base_idx = 0
		pos_idx = 2
		ref_idx = 1
		addX = -1.5

	# move base to origin
	pos = coord.copy()
	pos -= pos[base_idx]

	# rotate pos to be at pos x
	a_vector = pos[pos_idx] - pos[base_idx]
	b_vector = np.cross(a_vector, pos[ref_idx] - pos[base_idx])
	rot, _  = Rotation.match_vectors(np.array(((1., 0.,0.), (0., 1., 0.))), np.array((a_vector, b_vector)))
	pos = rot.apply(pos)
    
	# place Y
	addelem = 'H F Cl Br'.split()['ABCD'.index(ylabel)]
	try:
		os.mkdir('int-%s-%s' % (rxn, ylabel))
	except:
		continue
	with open('int-%s-%s/run.inp' % (rxn, ylabel)) as fh:
		fh.write('''! MP2 6-311G(d) COpt NoFrozenCore

%pal nprocs 1
end

%mp2 MaxCore 5000
end

%geom MaxIter 3000
end

*xyz -1 1
''')
		for idx in range(len(elem)):
			if idx == base_idx:
				fh.write('%s 0.0 $ 0.0 $ 0.0 $\n' % elem[idx])
				continue
			if idx == pos_idx:
				fh.write('%s %f $ 0.0 $ 0.0 $\n' % (pos[idx, 0], elem[idx]))
				continue
			fh.write('%s %f %f %f\n' % (elem[idx], *pos[idx]))
		fh.write('%s %f 0.0 $ 0.0 $\n*\n' % (addelem, addX))

elem, coord = read_xyz(sys.argv[1])
for RXN in 'e2 sn2'.split():
	for YLABEL in 'ABCD':
		do_it(elem, coord, RXN, YLABEL)

