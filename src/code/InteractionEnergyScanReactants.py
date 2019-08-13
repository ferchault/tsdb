import numpy as np
import pandas as pd
import re
import scipy.interpolate as sci
import sys

fn = sys.argv[1]

def clean_file(fn):
    """ Parse the summary file and create a pandas DF."""
    res = []
    regex = r"(?P<molecule>[^/]+)/gogp(?P<conformer>[^/]+)-../(?P<rxn>[^/-]+)-(?P<Y>[^/-]+)-(?P<dist>[^/-]+).*:.*  [ ]*(?P<energy>.*)"
    for line in open(fn):
        res.append(re.match(regex, line).groupdict())
    df = pd.DataFrame(res)
    df.dist = pd.to_numeric(df.dist, errors='coerce')
    df.conformer = pd.to_numeric(df.conformer, errors='coerce')
    df.energy = pd.to_numeric(df.energy, errors='coerce')
    return df
df = clean_file(fn)

def find_minima(inputdf):
    """ Takes explicit interaction energy scans and extracts spline minima"""
    results = []
    for name, group in df.groupby('molecule conformer rxn Y'.split()):
        molecule, conformer, rxn, Y = name
        s = group.sort_values(by='dist')
        try:
            spl = sci.UnivariateSpline(s.dist, s.energy, s=0, k=4)
        except:
            continue
        roots = spl.derivative().roots()
        if len(roots) != 1:
            continue
        root = roots[0]
        energy = float(spl(root))
        results.append({'Y': Y, 'conformer': conformer, 'molecule': molecule, 'rxn': rxn, 'dist': root, 'minimum': energy})
    return pd.DataFrame(results)
minima = find_minima(df)
minima['target'] = minima.apply(lambda row: row['molecule'].replace('0', row['Y']), axis=1)
df['target'] = df.apply(lambda row: row['molecule'].replace('0', row['Y']), axis=1)

# reactant E_D_B_E_A_B 3.8 MP2/6-311G(d)//MP2/6-311G(d) -632.163248264322
for idx, row in df.iterrows():
	print ('reactant scan %s %s %d %f MP2/6-311G(d)//MP2/6-311G(d) %f' % (row.rxn, row.target, row.conformer, row.dist, row.energy))

for ix, row in minima.iterrows():
	print ('reactant minimum %s %s %d %f MP2/6-311G(d)//MP2/6-311G(d) %f' % (row.rxn, row.target, row.conformer, row.dist, row.minimum))

