#!/usr/bin/env python
import sys
import glob
import cclib
from cclib.io import ccopen
import numpy as np
import logging
import enum
import os.path
import filecmp
import shutil
import numpy as np

# Input file templates
orca_ts_input = '''! MP2 6-311G(d) OptTS NoFrozenCore ExtremeSCF %s

%%pal nprocs 1
end

%%mp2 MaxCore 5000
end

*xyzfile %d %d inp.xyz
'''
orca_nf_input = '''! MP2 6-311G(d) NoFrozenCore ExtremeSCF NumFreq %s

%%pal nprocs 1
end

%%mp2 MaxCore 10000
end

%%geom
Calc_Hess true
NumHess true
end

*xyzfile %d %d inp.xyz
'''
orca_sd_input = '''! MP2 6-311G(d) NoFrozenCore ExtremeSCF Opt %s

%%geom
MaxIter 3000
inhess unit
MaxStep 0.01
Trust 0.1
end

%%mp2 MaxCore 5000
end

*xyzfile %d %d inp.xyz
'''

orca_go_input = '''! MP2 6-311G(d) Opt NoFrozenCore ExtremeSCF %s

%%pal nprocs 1
end

%%mp2 MaxCore 5000
end

*xyzfile %d %d inp.xyz
'''

def guess_bonds(atomtypes, coords):
	"""Based on MDAnalysis/topology/guessers.html#guess_bonds, which in turn uses the method of VMD. """
	import scipy.spatial.distance
	
	fudge_factor = 0.80
	vdwradii = {"C": 1.5, "F": 1.2, "H": 0.4, "N": 1.10, "O": 1.05, "S": 1.6, "P": 1.6, 'Cl': 1.0, 'Br': 1.15}
	lower_bound = 0.1
	max_vdw = max([vdwradii[t] for t in atomtypes])

	bonds = []
	for i, atom in enumerate(atomtypes):
		vdw_i = vdwradii[atomtypes[i]]
		max_d = (vdw_i + max_vdw) * fudge_factor

		dist = scipy.spatial.distance.cdist(coords[i:i+1], coords)[0]
		idx = np.where((dist > lower_bound) & (dist <= max_d))[0]

		for a in idx:
			j = a
			atom_j = atomtypes[j]

			if dist[a] < (vdw_i + vdwradii[atomtypes[j]]) * fudge_factor:
				if not (a, i) in bonds:
					bonds.append((i, a))
	return bonds

def is_single_molecule(bonds, numatoms):
	collected = [0]
	added = collected
	while len(added) != 0:
		added = []
		for atom in collected:
			for bond in bonds:
				if atom not in bond:
					continue
				other = (set(bond) - set((atom,))).pop()
				if other not in collected:
					added.append(other)
		
		collected += added
	return len(collected) == numatoms

class StageStates(enum.Enum):
	PENDING = 1
	INPROGRESS = 2
	COMPLETED = 4
	BROKEN = 5
	FAILED = 100

class Stage(object):
	def __init__(self, label, charge, spin, **kwargs):
		self._label, self._charge, self._spin = label, charge, spin
		self._options = kwargs
		self._get_folders_lastpath = None

	def get_label(self):
		return self._label

	def get_folders(self, path, label=None):
		if label is None:
			label = self.get_label()
		if path != self._get_folders_lastpath:
			self._get_folders_paths = [_ for _ in glob.glob('%s/*-*' % path) if os.path.isdir(_)]
		
		ret = sorted([_ for _ in self._get_folders_paths if _.startswith('%s/%s-' % (path.rstrip('/'), label))])
		return ret

	def _output_available(self, folder):
		return os.path.exists('%s/run.log' % folder)

	def run(self, path, req, verbose):
		self._verbose = verbose
		# Has anything been done, yet?
		folders = self.get_folders(path)

		if len(folders) == 0:
			fn = '%s/%s-00' % (path, self.get_label())
			return self._do_init(fn, path, req)

		# Is the last input fine?
		try:
			if not self._do_check_input(folders, req):
				self._verbose_message('Input check failed.')
				return StageStates.BROKEN
		except FileNotFoundError:
			self._verbose_message('No input file found.')
			return StageStates.BROKEN

		# Check output
		if self._output_available(folders[-1]):
			state = self._do_check_output(folders, req)
			if state is not None:
				return state
		else:
			return StageStates.INPROGRESS

		self._do_continue(folders, req)
		return StageStates.INPROGRESS

	def get_solvent(self):
		if self._options['implicit'] == 'gasphase':
			return ''
		if self._options['implicit'] == 'water':
			return 'CPCM(Water)'
		raise ValueError('Unknown solvent.')

	@staticmethod
	def quiet_ccopen(logfile):
		log = ccopen(logfile)
		logging.disable(logging.WARNING)
		data = log.parse()
		logging.disable(logging.NOTSET)
		return data

	@staticmethod
	def get_nextfolder(folders):
		parts = folders[-1].split('/')
		label, num = parts[-1].split('-')
		num = int(num) + 1
		return '%s/%s-%02d' % ('/'.join(parts[:-1]), label, num)

	def _verbose_message(self, message):
		if self._verbose:
			print ('%s: %s' % (self.get_label(), message))

	def _check_orca_input(self, lines, folders, cmds):
		try:
			solvent = self.get_solvent()
			try:
				extra = self._options['ignore_keywords']
			except:
				extra = []
			parts = cmds.split()
			if len(solvent) != 0:
				parts.append(solvent)
			actual = set(lines[0].strip().split()).union(set(extra))
			expected = set(parts)
			if actual != expected:
				self._verbose_message('Orca input line wrong.')
			assert(actual == expected)

			actual = lines[-1].strip()
			expected = '*xyzfile %d %d inp.xyz' % (self._charge, self._spin)
			if actual != expected:
				self._verbose_message('Spincharge line wrong.')
			assert (actual == expected)
		except AssertionError:
			return False

		if len(folders) > 1:
			lastrun, thisrun = folders[-2:]
			lastfile = '%s/run.xyz' % lastrun
			if not os.path.exists(lastfile):
				lastfile = '%s/inp.xyz' % lastrun
			try:
				if not filecmp.cmp(lastfile, '%s/inp.xyz' % thisrun):
					self._verbose_message('run.xyz from last run is not inp.xyz of this run')
					return False
			except FileNotFoundError:
				self._verbose_message('Previous run.xyz / current inp.xyz cannot be found.')
				return False

		return True

	def _xyz_continuation(self, folders, req):
		if len(req) != 1:
			self._verbose_message('No dependency defined, but continuation requested.')
		assert(len(req) == 1)
		path = '/'.join(folders[0].split('/')[:-1])
		reqfolders = self.get_folders(path, req[0])
		lastfile = '%s/run.xyz' % reqfolders[-1]
		if not os.path.exists(lastfile):
			lastfile = '%s/inp.xyz' % reqfolders[-1]
		if not filecmp.cmp(lastfile, '%s/inp.xyz' % folders[0]):
			self._verbose_message('Continuation: last run.xyz or current inp.xyz not found.')
			return False
		return True

	def _do_init_nf_scale(self, folder, path, req, scale):
		os.mkdir(folder)
		# Read NF data
		nflog = '%s/run.log' % (self.get_folders(path, req[0])[-1])
		data = self.quiet_ccopen(nflog)

		p = cclib.parser.utils.PeriodicTable()

		freqs = data.vibfreqs
		coords = data.atomcoords
		coords_TS = coords[-1]
		displacements = data.vibdisps
		
		atoms = [p.element[i] for i in data.atomnos]
		nAtoms = data.natom

		coords = coords_TS + displacements[0] * scale

		# Write new input files
		with open('%s/inp.xyz' % folder, 'w') as infile:
			infile.write('%d\n\n' % nAtoms)
			for i, atom in enumerate(coords):
				infile.write('%s %5.8f %5.8f %5.8f\n' % (atoms[i], *atom))

		with open('%s/run.inp' % folder, 'w') as runfile:
			runfile.write(orca_sd_input % (self.get_solvent(), self._charge, self._spin))

		return StageStates.INPROGRESS

class TS_Stage(Stage):
	def get_input(self):
		return orca_ts_input % (self.get_solvent(), self._charge, self._spin)

	def _do_init(self, folder, path, req):
		self._verbose_message('TS Stages cannot be initialised.')
		return StageStates.BROKEN

	def _do_check_input(self, folders, req):
		lines = open('%s/run.inp' % folders[-1], errors='replace').readlines()
		testA = self._check_orca_input(lines, folders, '! MP2 6-311G(d) OptTS NoFrozenCore ExtremeSCF')
		testB = self._check_orca_input(lines, folders, '! MP2 6-311G(d) OptTS NoFrozenCore')
		return testA or testB

	def _do_check_output(self, folders, req):
		loglines = open('%s/run.log' % folders[-1], errors='replace').readlines()

		# Validate log
		orca_ok, maxstep, geoopt_complete = False, False, False
		for line in loglines:
			if 'ORCA TERMINATED NORMALLY' in line: orca_ok = True
			if 'The optimization did not converge but reached the maximum number of' in line: maxstep = True
			if 'THE OPTIMIZATION HAS CONVERGED' in line: geoopt_complete = True

		if orca_ok and geoopt_complete and not maxstep:
			return StageStates.COMPLETED

	def _do_continue(self, folders, req):
		nextfolder = Stage.get_nextfolder(folders)
		os.mkdir(nextfolder)

		lastfolder = folders[-1]
		shutil.copyfile('%s/run.xyz' % lastfolder, '%s/inp.xyz' % nextfolder)
		shutil.copyfile('%s/run.inp' % lastfolder, '%s/run.inp' % nextfolder)

class NF_Stage(Stage):
	@staticmethod
	def validate_ignore(logfile):
		return True

	@staticmethod
	def validate_e2_ts(logfile):
		def check_neg_freqs(freqs):
			return np.amin(freqs) < -400 and len(freqs[freqs < 0]) == 1

		def check_displacements_e2(coords, displacements):
			# get bonds
			bond_LG = coords[0] - coords[2]
			bond_LG_new = (coords[0] + displacements[0]) - (coords[2] + displacements[2])
			bond_H = coords[1] - coords[4]
			bond_H_new = (coords[1] + displacements[1]) - (coords[4] + displacements[4])
			bond_Nu = coords[4] - coords[3]
			bond_Nu_new = (coords[4] + displacements[4]) - (coords[3] + displacements[3])

			# check displacements 
			if np.linalg.norm(coords[1] - coords[4]) < 1.15:
				return False
			if (np.linalg.norm(bond_LG) > np.linalg.norm(bond_LG_new)):
				LG = 1
			else:
				LG = 0
			if (np.linalg.norm(bond_H) > np.linalg.norm(bond_H_new)):
				H = 1
			else:
				H = 0
			if (np.linalg.norm(bond_Nu) > np.linalg.norm(bond_Nu_new)):
				Nu = 1
			else:
				Nu = 0

			return (LG == 1 and H == 1 and Nu == 0) or (LG == 0 and H == 0 and Nu == 1)

		def check_CX_vs_XH(labels, coords):
			C1 = coords[0]
			X  = coords[2]
			CX = np.array([ C1[0] - X[0], C1[1] - X[1], C1[2] - X[2] ])
			CX_len = np.linalg.norm(CX)

			for i in range(len(labels)):
				if labels[i] == 1:
					H = coords[i]
					XH = np.array([ H[0] - X[0], H[1] - X[1], H[2] - X[2] ])
					XH_len = np.linalg.norm(XH)

					if CX_len > XH_len: return False
			return True

		data = Stage.quiet_ccopen(logfile)
		
		return check_neg_freqs(data.vibfreqs) and check_displacements_e2(data.atomcoords[-1], data.vibdisps[0]) and check_CX_vs_XH(data.atomnos, data.atomcoords[-1])

	@staticmethod
	def validate_sn2_ts(logfile):
		def check_neg_freqs(freqs):
			return np.amin(freqs) < -400 and len(freqs[freqs < 0]) == 1

		def check_displacements_sn2(coords, displacements):
			# get bonds
			bond_LG = coords[0] - coords[2]
			bond_LG_new = (coords[0] + displacements[0]) - (coords[2] + displacements[2])
			bond_Nu = coords[0] - coords[3]
			bond_Nu_new = (coords[0] + displacements[0]) - (coords[3] + displacements[3])

			# check displacements 
			if (np.linalg.norm(bond_LG) > np.linalg.norm(bond_LG_new)):
				LG = 1
			else:
				LG = 0
			if (np.linalg.norm(bond_Nu) > np.linalg.norm(bond_Nu_new)):
				Nu = 1
			else:
				Nu = 0

			return (LG == 1 and Nu == 0) or (LG == 0 and Nu == 1)

		data = Stage.quiet_ccopen(logfile)

		return check_neg_freqs(data.vibfreqs) and check_displacements_sn2(data.atomcoords[-1], data.vibdisps[0])

	def _do_init(self, folder, path, req):
		os.mkdir(folder)
		with open('%s/run.inp' % folder, 'w') as runfile:
			runfile.write(orca_nf_input % (self.get_solvent(), self._charge, self._spin))

		assert (len(req) == 1)
		reqfolders = self.get_folders(path, req[0])
		shutil.copyfile('%s/run.xyz' % reqfolders[-1], '%s/inp.xyz' % folder)

		return StageStates.INPROGRESS

	def _do_check_input(self, folders, req):
		lines = open('%s/run.inp' % folders[-1], errors='replace').readlines()
		testA = self._check_orca_input(lines, folders, '! MP2 6-311G(d) OptTS NoFrozenCore ExtremeSCF NumFreq')
		testB = self._check_orca_input(lines, folders, '! MP2 6-311G(d) NoFrozenCore ExtremeSCF NumFreq')
		testC = self._check_orca_input(lines, folders, '! MP2 6-311G(d) NoFrozenCore NumFreq')
		if not testA and not testB and not testC:
			return False

		if self._options['continuation_check']:
			if not self._xyz_continuation(folders, req):
				return False

		return True

	def _do_continue(self, folders, req):
		nextfolder = Stage.get_nextfolder(folders)
		os.mkdir(nextfolder)

		lastfolder = folders[-1]
		shutil.copyfile('%s/inp.xyz' % lastfolder, '%s/inp.xyz' % nextfolder)
		shutil.copyfile('%s/run.inp' % lastfolder, '%s/run.inp' % nextfolder)

		foundhess = False
		for folder in folders[::-1]:
			for fn in glob.glob('%s/*.hess' % lastfolder):
				foundhess = True
				shutil.copyfile(fn, '%s/%s' % (nextfolder, os.path.basename(fn)))
			if foundhess:
				break

	def _do_check_output(self, folders, req):
		fn = '%s/run.log' % folders[-1]
		loglines = open(fn, errors='replace').readlines()

		# Validate log
		orca_ok = False
		for line in loglines:
			if 'ORCA TERMINATED NORMALLY' in line: orca_ok = True

		# Check presence of .hess files
		if not orca_ok:
			numfiles = len(glob.glob('%s/*.hess' % folders[-1]))
			if numfiles == 0:
				return StageStates.FAILED

		# Validate TS
		if orca_ok:
			orca_ts_ok = self._options['validation'](fn)

		if orca_ok and orca_ts_ok:
			return StageStates.COMPLETED

		if orca_ok and not orca_ts_ok:
			return StageStates.FAILED

class LS_Stage(Stage):
	def _do_check_input(self, folders, req):
		lines = open('%s/run.inp' % folders[-1], errors='replace').readlines()
		if not self._check_orca_input(lines, folders, '! MP2 6-311G(d) Opt NoFrozenCore ExtremeSCF'):
			return False
		return True

	def _do_check_output(self, folders, req):
		orca_ok = False
		for line in open('%s/run.log' % folders[-1], errors='replace').readlines():
			if 'ORCA TERMINATED NORMALLY' in line: break
		else:
			return StageStates.FAILED
		return StageStates.COMPLETED

	def _do_init(self, folder, path, req):
		if self._options['forward']:
			scale = 0.01 * self._options['step']
		else:
			scale = -0.01 * self._options['step']
		return self._do_init_nf_scale(folder, path, req, scale)

class SD_Stage(Stage):
	@staticmethod
	def validate_single_molecule(xyzfile):
		lines = open(xyzfile).readlines()[2:]
		elements = [_.strip().split()[0] for _ in lines]
		coords = np.array([list(map(float, _.strip().split()[1:])) for _ in lines])
		bonds = guess_bonds(elements, coords)
		return is_single_molecule(bonds, len(elements))

	def _do_check_input(self, folders, req):
		lines = open('%s/run.inp' % folders[-1], errors='replace').readlines()
		if not self._check_orca_input(lines, folders, '! MP2 6-311G(d) Opt NoFrozenCore ExtremeSCF'):
			return False
		return True

	def _do_check_output(self, folders, req):
		orca_ok, geoopt_complete = False, False
		for line in open('%s/run.log' % folders[-1], errors='replace').readlines():
			if 'ORCA TERMINATED NORMALLY' in line: orca_ok = True
			if 'THE OPTIMIZATION HAS CONVERGED' in line: geoopt_complete = True
			if orca_ok and geoopt_complete: break
		else:
			if not os.path.exists('%s/run.xyz' % folders[-1]):
				return StageStates.FAILED
			return None
		if 'validation' in self._options:
			xyzfile = '%s/run.xyz' % folders[-1]
			if not self._options['validation'](xyzfile):
				return StageStates.FAILED
		return StageStates.COMPLETED

	def _do_continue(self, folders, req):
		nextfolder = Stage.get_nextfolder(folders)
		os.mkdir(nextfolder)

		lastfolder = folders[-1]
		shutil.copyfile('%s/run.xyz' % lastfolder, '%s/inp.xyz' % nextfolder)
		shutil.copyfile('%s/run.inp' % lastfolder, '%s/run.inp' % nextfolder)

	def _do_init_nf(self, folder, path, req):
		if self._options['forward']:
			scale = 0.1
		else:
			scale = -0.1
		return self._do_init_nf_scale(folder, path, req, scale)

	def _do_init_sp(self, folder, path, req):
		filenames = {'R': 'reactant.xyz', 'P': 'product.xyz'}
		index = self.get_label()[-1]

		os.mkdir(folder)
		reqfolder = self.get_folders(path, req[0])[-1]

		shutil.copyfile('%s/%s' % (reqfolder, filenames[index]), '%s/inp.xyz' % folder)
		with open('%s/run.inp' % folder, 'w') as runfile:
			runfile.write(orca_go_input % (self.get_solvent(), self._charge, self._spin))

		return StageStates.INPROGRESS

	def _do_init(self, folder, path, req):
		if req[0].startswith('nf'):
			return self._do_init_nf(folder, path, req)
		if req[0].startswith('sp'):
			return self._do_init_sp(folder, path, req)
		raise NotImplementedError('No such parent known.')

class Split_Stage(Stage):
	@staticmethod
	def _identify_type(folder, rxn):
		p = cclib.parser.utils.PeriodicTable()

		data = Stage.quiet_ccopen('%s/run.log' % folder)

		coords, nAtoms = data.atomcoords, data.natom
		atoms = [p.element[i] for i in data.atomnos]

		# take last geometry
		coords = coords[-1]

		# check bond lengths
		if rxn == 'sn2':
			if np.linalg.norm(coords[0] - coords[3]) >  np.linalg.norm(coords[0] - coords[2]):
				coords = np.delete(coords, 3, 0)
				atoms = np.delete(atoms, 3, 0)
				side = "reactant"
			elif np.linalg.norm(coords[0] - coords[3]) <  np.linalg.norm(coords[0] - coords[2]):
				coords = np.delete(coords, 2, 0)
				atoms = np.delete(atoms, 2, 0)
				side = "product"

		if rxn == 'e2':
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

		xyz = ''
		if side == "product":
			xyz += '%s\n\n' % (nAtoms - 3)
		else:
			xyz += '%s\n\n' % (nAtoms - 1)

		for i, atom in enumerate(coords):
			xyz += '%s %5.8f %5.8f %5.8f\n' % (str(atoms[i]), *atom)

		return side, xyz

	def _output_available(self, folder):
		return True

	def _do_init(self, folder, path, req):
		os.mkdir(folder)
		assert (len(req) == 2)
		reqfolders = [self.get_folders(path, _) for _ in req]
		identified = [Split_Stage._identify_type(_[-1], self._options['rxn']) for _ in reqfolders]
		
		# Exactly one R and one P?
		types = [_[0] for _ in identified]
		if set(types) != set('reactant product'.split()):
			self._verbose_message('Not both reactant and product found.')
			return StageStates.BROKEN

		# Write new structures to xyz files
		for moltype, xyz in identified:
			with open('%s/%s.xyz' % (folder, moltype), 'w') as r:
				r.write(xyz)

		return StageStates.COMPLETED

	def _do_check_input(self, folders, req):
		return True

	def _do_check_output(self, folders, req):
		for moltype in 'reactant product'.split():
			if not os.path.exists('%s/%s.xyz' % (folders[-1], moltype)):
				self._verbose_message('Not both reactant and product found.')
				return StageStates.BROKEN
		return StageStates.COMPLETED

class Workflow(object):
	def __init__(self):
		self._stages = []
		self._labels = []
		self._requirements = []

	def add_stage(self, stage, required):
		label = stage.get_label()
		if label in self._labels:
			raise ValueError('Label not unique for this workflow.')
		self._stages.append(stage)
		self._labels.append(label)
		self._requirements.append(required)

	def show_dependencies(self):
		print ('digraph G {')
		for reqs, stage in zip(self._requirements, self._stages):
			for req in reqs:
				print ('"%s" -> "%s"' % (req, stage.get_label()))
		print ('}')

	def read_path(self, path, verbose=False):
		progress = {state: [] for state in StageStates}
		while True:
			runnable = False
			for req, stage in zip(self._requirements, self._stages):
				if stage.get_label() in sum(progress.values(), []):
					continue

				if set(req).issubset(set(progress[StageStates.COMPLETED])):
					runnable = True
					state = stage.run(path, req, verbose=verbose)
					progress[state].append(stage.get_label())
					break
			if not runnable:
				break

		progress[StageStates.PENDING] = list(set(self._labels) - set(sum(progress.values(), [])))
		if verbose:
			for state in StageStates:
				print (state, ' '.join(sorted(progress[state])))
		return len(progress[StageStates.COMPLETED]) == len(self._stages), progress

class E2Workflow(Workflow):
	def __init__(self, implicit='gasphase'):
		super().__init__()
		self.add_stage(TS_Stage(label='ts', charge=-1, spin=1, implicit=implicit), [])
		self.add_stage(NF_Stage(label='nf', charge=-1, spin=1, implicit=implicit, validation=NF_Stage.validate_e2_ts, continuation_check=False), ['ts'])
		for step in range(1, 11):
			self.add_stage(LS_Stage(label='lsA%s' % (step - 1), charge=-1, spin=1, implicit=implicit, forward=True, step=step), ['nf'])
			self.add_stage(LS_Stage(label='lsB%s' % (step - 1), charge=-1, spin=1, implicit=implicit, forward=False, step=step), ['nf'])
		self.add_stage(SD_Stage(label='sdA', charge=-1, spin=1, implicit=implicit, forward=True), ['nf'])
		self.add_stage(SD_Stage(label='sdB', charge=-1, spin=1, implicit=implicit, forward=False), ['nf'])
		self.add_stage(Split_Stage(label='sp', charge=-1, spin=1, implicit=implicit, rxn='e2'), ['sdA', 'sdB'])
		for label in 'RP':
			self.add_stage(SD_Stage(label='go%s' % label, charge=0, spin=1, implicit=implicit, validation=SD_Stage.validate_single_molecule), ['sp'])
			self.add_stage(NF_Stage(label='nf%s' % label, charge=0, spin=1, implicit=implicit, continuation_check=True, validation=NF_Stage.validate_ignore), ['go%s' % label])

class E2LightWorkflow(Workflow):
	def __init__(self, implicit='gasphase'):
		super().__init__()
		self.add_stage(TS_Stage(label='ts', charge=-1, spin=1, implicit=implicit, ignore_keywords=['ExtremeSCF']), [])
		self.add_stage(NF_Stage(label='nf', charge=-1, spin=1, implicit=implicit, validation=NF_Stage.validate_e2_ts, continuation_check=False), ['ts'])
		self.add_stage(SD_Stage(label='sdA', charge=-1, spin=1, implicit=implicit, forward=True), ['nf'])
		self.add_stage(SD_Stage(label='sdB', charge=-1, spin=1, implicit=implicit, forward=False), ['nf'])
		self.add_stage(Split_Stage(label='sp', charge=-1, spin=1, implicit=implicit, rxn='e2'), ['sdA', 'sdB'])

class SN2Workflow(Workflow):
	def __init__(self, implicit='gasphase'):
		super().__init__()
		self.add_stage(TS_Stage(label='ts', charge=-1, spin=1, implicit=implicit), [])
		self.add_stage(NF_Stage(label='nf', charge=-1, spin=1, implicit=implicit, validation=NF_Stage.validate_sn2_ts, continuation_check=False), ['ts'])
		for step in range(1, 11):
			self.add_stage(LS_Stage(label='lsA%s' % (step - 1), charge=-1, spin=1, implicit=implicit, forward=True, step=step), ['nf'])
			self.add_stage(LS_Stage(label='lsB%s' % (step - 1), charge=-1, spin=1, implicit=implicit, forward=False, step=step), ['nf'])
		self.add_stage(SD_Stage(label='sdA', charge=-1, spin=1, implicit=implicit, forward=True), ['nf'])
		self.add_stage(SD_Stage(label='sdB', charge=-1, spin=1, implicit=implicit, forward=False), ['nf'])
		self.add_stage(Split_Stage(label='sp', charge=-1, spin=1, implicit=implicit, rxn='sn2'), ['sdA', 'sdB'])
		for label in 'RP':
			self.add_stage(SD_Stage(label='go%s' % label, charge=0, spin=1, implicit=implicit, validation=SD_Stage.validate_single_molecule), ['sp'])
			self.add_stage(NF_Stage(label='nf%s' % label, charge=0, spin=1, implicit=implicit, continuation_check=True, validation=NF_Stage.validate_ignore), ['go%s' % label])


if __name__ == '__main__':
	rxn, solvent, path = sys.argv[1:]
	if rxn == 'sn2':
		wf = SN2Workflow(solvent)
	if rxn == 'e2':
		wf = E2Workflow(solvent)
	#import cProfile, pstats, io
	#pr = cProfile.Profile()
	#pr.enable()
	wf.read_path(path, verbose=True)
	#pr.disable()
	#s = io.StringIO()
	#ps = pstats.Stats(pr, stream=s)
	#ps.sort_stats('cumulative')
	#ps.print_stats()
	#print(s.getvalue())
		
