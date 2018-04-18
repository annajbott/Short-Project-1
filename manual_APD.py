#!/usr/bin/env python

import matplotlib.pyplot as pl
import myokit
import numpy as np

#Set pacing
pcl = 200

# Load model and set protocol, create simulation
m = myokit.load_model('ohara-cipa-v1-2017.mmt')
p = myokit.pacing.blocktrain(pcl, 0.5, offset=0, level=1.0, limit=0)
s = myokit.Simulation(m,p)

# Set cell type
cell_types = {'Endocardial' : 0, 'Epicardial' : 1, 'Mid-myocardial' : 2}
cell_type =  'Epicardial'
s.set_constant('cell.celltype', cell_types[cell_type])

# Pre-pace simulation to steady state
s.pre(400*pcl)
d = s.run(pcl)
potential = np.asarray(d['membrane.V'])
min_potential = np.ndarray.min(potential)
max_potential = np.ndarray.max(potential)
range_AP = max_potential - min_potential
range10 = 0.1*range_AP
thresh = min_potential + range10
print min_potential, max_potential, thresh
s.reset()
paces = 5
d = s.run(paces*pcl)
potential = np.asarray(d['membrane.V'])
time = np.asarray(d['environment.time'])

# Assuming no 1:2 (or greater) ratios, so will not be more starts than number of paces
start_apd = np.zeros(paces)
duration_apd = np.zeros(paces)
status = np.zeros(len(potential))

# Recording start and duration of APDs
in_apd = 0
responses = 0
for i in range(0,len(potential)):
    # Start of AP, greater than threshold
    if potential[i] >= thresh and not in_apd:
        in_apd = 1
        start_apd[responses] = time[i]
        status[i] = 1
    # In actional potential, in_apd = 1, still greater than threshold
    elif in_apd and potential[i] >= thresh:
        status[i] = 1

    # End of AP, record end time and subtract start for the duration
    elif in_apd and potential[i] < thresh:
        in_apd = 0
        end_apd = time[i]
        duration_apd[responses] = end_apd - start_apd[responses]
        status[i] = 0
        responses += 1
    else:
        status[i] = 0

print start_apd
print duration_apd








# Plot the results
pl.figure()
pl.plot(d['environment.time'], d['membrane.V'])
pl.show()
