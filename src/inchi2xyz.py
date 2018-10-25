#!/usr/bin/env python
from rdkit import Chem
from rdkit.Chem import AllChem

import sys
mp = Chem.AddHs(Chem.MolFromInchi(sys.argv[1]))
AllChem.EmbedMolecule(mp)
AllChem.UFFOptimizeMolecule(mp)
coords = mp.GetConformer().GetPositions()
elements = [_.GetSymbol() for _ in mp.GetAtoms()]
print (len(coords))
print ('')
for label, pos in zip(elements, coords):
	print (label, pos[0], pos[1], pos[2])
