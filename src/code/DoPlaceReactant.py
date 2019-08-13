#!/usr/bin/env python
import sys
import os
import numpy as np

lines = [_.strip() for _ in open(sys.argv[1] + '/run.xyz').readlines()[2:]]
elements = [_.split()[0] for _ in lines]
coords = np.array([list(map(float, _.split()[1:])) for _ in lines])
Y = {'A': 'H', 'B': 'F', 'C': 'Cl', 'D': 'Br'}[sys.argv[2]]
ds = [float(sys.argv[3])]
rxn = sys.argv[4]

for d in ds:
	# calculate position
	c1 = coords[0]
	c2 = coords[1]
	X = coords[2]
	H = coords[3]

	# e2
	tcoord = dict()
	direction = H - c2
	direction /= np.linalg.norm(direction)
	tcoord['e2'] = c2 + d * direction

	# sn2
	direction = c1 - X
	direction /= np.linalg.norm(direction)
	tcoord['sn2'] = c1 + d * direction

	coord = tcoord[rxn]

	# write input geometry
	print ('%d\n' % (len(lines) + 1))
	for line in lines:
		print ('%s' % line)
	print ('%s %s %s %s' % (Y, coord[0], coord[1], coord[2]))

