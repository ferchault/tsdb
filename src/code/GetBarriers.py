#!/usr/bin/env python

import pandas as pd
import numpy as np
import scipy.interpolate as sci

basedir='RESULTS/'

# Get valid calculations (i.e. have the sp stage completed for that calc)
completed = pd.read_csv(basedir + 'completed.txt', sep=' ', header=None, names='rxn calc step'.split())
completed = completed.query('step=="nf"')

# Obtain TS energies for reference level of theory
tsenergies = pd.read_csv(basedir + 'total-electronic-energies.txt', sep=' ', header=None, names='rxn calc target lot step energy'.split())
tsenergies = tsenergies.query('step=="ts" & lot=="MP2/6-311G(d)//MP2/6-311G(d)"')

# Filter for validated TS
tsenergies = pd.merge(completed, tsenergies, on='rxn calc'.split(), how='left')
tsenergies = tsenergies['rxn calc target energy'.split()]

# Select lowest TS for each target
tsenergies = tsenergies.sort_values('energy')
tsenergies['ts_rank'] = pd.to_numeric(tsenergies.groupby('target rxn'.split())['energy'].rank(method='dense'), downcast='integer')

# Get interaction minima
interaction = pd.read_csv(basedir + 'interaction-energies.txt', sep=' ', header=None, names='geotype method rxn target conf dist lot energy'.split())
interaction = interaction.query('method=="minimum" & lot=="MP2/6-311G(d)//MP2/6-311G(d)"')
interaction = interaction.sort_values('energy')
interaction['conformer_rank'] = pd.to_numeric(interaction.groupby('target rxn'.split())['energy'].rank(method='dense'), downcast='integer')
interaction = interaction['rxn target conf dist lot energy conformer_rank'.split()]

# Calculate barriers
barriers = pd.merge(tsenergies, interaction, on='rxn target'.split(), how='inner', suffixes='_ts _r'.split())
barriers['barrier'] = barriers.energy_ts - barriers.energy_r

# Write output
barriers.to_csv(basedir + 'electronic-forward-barriers.txt', sep=' ', header=False, index=False)

