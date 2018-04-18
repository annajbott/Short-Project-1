#!/usr/bin/env python

import matplotlib.pyplot as pl
import myokit
import numpy as np

#Set pacing
pcl = 150

# Load model and set protocol, create simulation
m = myokit.load_model('ohara-cipa-v1-2017.mmt')
p = myokit.pacing.blocktrain(pcl, 0.5, offset=30, level=1.0, limit=0)
s = myokit.Simulation(m,p)

# Set cell type
cell_types = {'Endocardial' : 0, 'Epicardial' : 1, 'Mid-myocardial' : 2}
cell_type =  'Epicardial'
s.set_constant('cell.celltype', cell_types[cell_type])

# Pre-pace simulation to steady state and run 5 cycles of simulation
s.pre(400*pcl)
d = s.run(5*pcl)

# Convert membrane potential arrays to numpy
potential = np.asarray(d['membrane.V'])

# Set threshold, 10% of range of min max membrane potentials
min_potential = np.ndarray.min(potential)
max_potential = np.ndarray.max(potential)
range_10AP = 0.1*(max_potential- min_potential)
thresh = min_potential + range_10AP

# Check min and max are reasonable
print min_potential, max_potential, thresh

# Reset simulation after 2 runs to obtain threshold
s.reset()

# Run actual simulation to calculate APD for
paces = 20
d = s.run(paces*pcl )
potential = np.asarray(d['membrane.V'])
time = np.asarray(d['environment.time'])

# Numpy arrays containing start of APs and duration of APs
# Assuming no 1:2 (or 1:3 etc) ratios, so number of paces >= number of APs
start_ap = np.zeros(paces)
duration_ap = np.zeros(paces)

in_apd = 0 # Initial is to not be in an AP
responses = 0
# Recording start and duration of APDs
# Assuming end of AP is when it is below threshold (~75mV), could measure max for
for i in range(0,len(potential)):
    # Start of AP, greater than threshold
    if potential[i] >= thresh and not in_apd:
        in_apd = 1
        start_ap[responses] = time[i]
    # In actional potential, in_apd = 1, still greater than threshold
    #elif in_apd and potential[i] >= thresh:


    # End of AP, record end time and subtract start for the duration
    elif in_apd and potential[i] < thresh:
        in_apd = 0
        end_apd = time[i]
        duration_ap[responses] = end_apd - start_ap[responses]

        responses += 1

# If started simulation in the middle of an AP, get rid of first entry
if potential[0] > thresh:
    start_ap = start_ap[1:]
    duration_ap = duration_ap[1:]
print start_ap
print duration_ap








# Plot the results
pl.figure()
pl.plot(d['environment.time'], d['membrane.V'])
pl.show()
