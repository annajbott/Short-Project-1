#!/usr/bin/env python

import matplotlib.pyplot as pl
import myokit
import numpy as np

def ap_duration(s, bcl, paces):
    d = s.run(paces*bcl )
    # Convert membrane potential and time arrays to numpy
    potential = np.asarray(d['membrane.V'])
    time = np.asarray(d['environment.time'])

    # Numpy arrays containing start of APs and duration of APs
    # Assuming no 1:2 (or 1:3 etc) ratios, so number of paces >= number of APs
    start_ap = np.zeros(paces + 1) # Plus one in case starts mid AP
    duration_ap = np.zeros(paces + 1)

    # Recording start and duration of APDs
    # Assuming end of AP is when it is below threshold (~75mV), could measure max and calculate 90% of each AP
    # Initialising. Assuming not in AP, will delete first entry if starts in AP
    in_ap = 0
    responses = 0
    for i in range(0,len(potential)):
        # Start of AP (so not already in AP), greater than threshold
        if potential[i] >= thresh and not in_ap:
            in_ap = 1
            start_ap[responses] = time[i]
        # End of AP, record end time and subtract start for the duration
        elif in_ap and potential[i] < thresh:
            in_ap = 0
            end_apd = time[i]
            duration_ap[responses] = end_apd - start_ap[responses]
            responses += 1

    # If simulation was started in the middle of an AP, discard first entry
    if potential[0] > thresh:
        start_ap = start_ap[1:]
        duration_ap = duration_ap[1:]
    # If final AP not finished, discard its start time
    if potential[-1] > thresh:
        start_ap = start_ap[0:-1]

    # Keep only non-zero terms
    start_ap = start_ap[np.nonzero(start_ap)]
    duration_ap = duration_ap[np.nonzero(duration_ap)]

    # Plot the results
    pl.figure()
    pl.plot(d['environment.time'], d['membrane.V'])

    return[start_ap, duration_ap]

#Set pacing
bcl = 50

# Load model and set protocol, create simulation
m = myokit.load_model('ohara-cipa-v1-2017.mmt')
p = myokit.pacing.blocktrain(bcl, 0.5, offset=0, level=1.0, limit=0)
s = myokit.Simulation(m,p)

# Set cell type
cell_types = {'Endocardial' : 0, 'Epicardial' : 1, 'Mid-myocardial' : 2}
cell_type =  'Epicardial'
s.set_constant('cell.celltype', cell_types[cell_type])

# Pre-pace simulation to steady state and run 5 cycles of simulation
s.pre(40*bcl)

# Set threshold, 90% of repolarisation to resting membrane potential
thresh = 0.9 * s.state()[m.get('membrane.V').indice()]

# Run actual simulation to calculate APD for
paces = 20

# Running using function
start, duration = ap_duration(s,bcl, paces)
print start, duration
pl.show()
