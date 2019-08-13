import cclib
import glob
import warnings
warnings.filterwarnings("ignore")

def read_fn(prefix, fn):
	try:
		data = cclib.io.ccread(fn)
		for kind in data.atomcharges:
			for idx in range(len(data.atomcharges[kind])):
				print (prefix, kind, idx, data.atomcharges[kind][idx])
	except:
		pass

# reactants
def do_reactants():
	for folder in glob.glob('/mnt/ECHO/bookie/conformer/reactants/*/'):
		label = folder.strip('/').split('/')[-1]
		for fn in glob.glob('%s/gogp*-00/run.log' % folder):
			i = int(fn.split('/')[-2].split('-')[0][-2:])
			read_fn('reactant %s %s %s' % (label, label, i), fn)
do_reactants()

def clean_target(target):
	target = target.replace('e2ts_', '')
	target = target.replace('sn2_', '')
	if '-' not in target and len(target) > len('A_A_A_A_A_A'):
		target = target[:len('A_A_A_A_A_A')] + '-' + target[-len('A_A_A_A_A_A'):]
	return target


def do_ts(rxn):
	for folder in glob.glob('/mnt/ECHO/bookie/%s/*/' % rxn):
		try:
			calc = clean_target(folder.strip('/').split('/')[-1])
			target = calc[-len('A_A_A_A_A_A'):]
			steps = sorted(glob.glob('%s/ts-??/' % folder))
			if len(steps) == 0:
				continue
			step = steps[-1]
			fn = '%s/run.log' % (step)
			read_fn('ts %s %s %s' % (rxn, calc, target), fn)
		except:
			continue

do_ts('e2')
do_ts('sn2')

