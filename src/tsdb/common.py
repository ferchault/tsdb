#!/usr/bin/env python
# -*- coding: utf-8 -*-
substituents = {'A': 'H'.split(), 'B': 'N O O'.split(), 'C': 'C N'.split(), 'D': 'C H H H'.split(), 'E': 'N H H'.split()}
leaving_X = {'A': 'F', 'B': 'Cl', 'C': 'Br'}
attack_Y = {'A': 'H', 'B': 'F', 'C': 'Cl', 'D': 'Br'}

def get_elements_from_label(label):
	""" Returns a set of elements from the label of the configuration alone."""
	parts = label.split('_')

	if len(parts) != 7:
		raise ValueError('Ill-formatted label')

	rxn = parts[0]
	if rxn not in 'e2ts sn2'.split():
		raise ValueError('Unkown rxn.')
	X = parts[5]
	Y = parts[6]

	if rxn == 'sn2':
		elements = 'C H'.split()

		elements += [leaving_X[X]]
		elements += [attack_Y[Y]]
		elements += ['C']

	if rxn == 'e2ts':
		elements = 'C C'.split()

		elements += [leaving_X[X]]
		elements += [attack_Y[Y]]
		elements += ['H']

	for R in parts[1:5]:
		elements += substituents[R]

	return sorted(elements)

def get_elements_from_xyz(lines):
	""" Obtains the set of elements from the lines of a xyz file."""
	elements = [_.strip().split()[0] for _ in lines[2:]]
	return sorted(elements)
