#!/usr/bin/env python

import matplotlib.pyplot as pl
import myokit
import numpy as np

def ap_duration(d, paces, threshold):

    # Convert membrane potential and time arrays to numpy
    potential = np.asarray(d['membrane.V'])
    time = np.asarray(d['environment.time'])

    # Numpy arrays containing start of APs and duration of APs
    # Assuming no 1:2 (or 1:3 etc) ratios, so number of paces >= number of APs
    start_ap = np.zeros(paces + 2) # Plus 2 in case starts mid AP and ends with incomplete AP
    duration_ap = np.zeros(paces + 2)

    # Recording start and duration of APDs
    # Assuming end of AP is when it is below threshold (~75mV), could measure max and calculate 90% of each AP
    # Initialising. Assuming not in AP, will delete first entry if starts in AP
    in_ap = 0
    responses = 0
    for i in range(0,len(potential)):
        # Start of AP (so not already in AP), greater than threshold
        if potential[i] >= threshold and not in_ap:
            in_ap = 1
            start_ap[responses] = time[i]
        # End of AP, record end time and subtract start for the duration
        elif in_ap and potential[i] < threshold:
            in_ap = 0
            end_apd = time[i]
            duration_ap[responses] = end_apd - start_ap[responses]
            responses += 1

    # If simulation was started in the middle of an AP, discard first entry
    if potential[0] > threshold:
        start_ap = start_ap[1:]
        duration_ap = duration_ap[1:]
    # If final AP not finished, discard its start time
    if potential[-1] > threshold:
        start_ap = start_ap[0:-1]

    # Keep only non-zero terms
    start_ap = start_ap[np.nonzero(start_ap)]
    duration_ap = duration_ap[np.nonzero(duration_ap)]

    # Plot the results
    #pl.figure()
    #pl.plot(d['environment.time'], d['membrane.V'])

    return[start_ap, duration_ap]

def steady(m,p,pcl,paces):
    # Create simulation using user defined model and protocol
    s = myokit.Simulation(m,p)

    # Progress report
    print "In steady function"

    # Run for pcl*paces, specifying to use time values every 0.1 ms interval for solver
    d  = s.run(paces*pcl, log_interval = 0.1)

    # Lengthy step due to log_interval component, so progress report when this is finished
    print "Simulations finished running"

    # Set repolarisation threshold to 90% and run APD function to get information
    # about when each AP starts
    thresh = 0.9 * s.state()[m.get('membrane.V').indice()]
    start, duration = ap_duration(d, paces, thresh)

    # If there is greater than a 1:1 ratio
    paces_per_period = np.round(float(start[-1] - start[-2]),0)/ float(pcl)
    # How many time points in each period (AP) (time interval for solver 0.1)
    time_points_per_period = pcl*paces_per_period/0.1

    # Convert membrane potential list into numpy array
    V = np.asarray(d['membrane.V'])

    # If multiple paces per period, V may not be divisible by no. time_points_per_period
    # Slice the end of the array off so it is divisble and can reshape array
    remainder = int(len(V)%time_points_per_period)

    if remainder != 0:
        V = V[0:-remainder]

    # Reshape array so that each period is on a seperate row
    # Each column contains the same point in an AP, but for different beats
    V = np.reshape(V, (-1, int(time_points_per_period)))

    # Euclidean distance and percentage change between adjacent periods calculated
    perc = np.zeros(len(V)-1)
    # Length V = number of rows --> number of APs (periods)
    for i in range(1,len(V)):
        dist = np.sqrt(np.sum((V[i]-V[i-1])**2))
        perc[i-1] = (dist/abs(np.mean(V[i-1])))*100

    # Moving average for 10 terms calculated for percentage change
    moving_av = np.convolve(perc, np.ones((10,))/10, mode='valid')

    # Plot results
    pl.figure()
    pl.plot(range(1,len(V)), perc)
    pl.xlabel('Period Number')
    pl.ylabel('Euclidean distance percentage change to previous period')
    pl.figure()
    pl.plot(range(1,len(V)-9), moving_av)
    pl.xlabel('Period Number')
    pl.ylabel('Moving average')
    print moving_av
    if np.size(np.nonzero(moving_av < 1)) < 1:
        ss = "Run again with more paces"
    else:
        ss = paces_per_period*np.nonzero(moving_av < 1)[0][0]
    return(ss)

def main():
    #Set pacing
    bcl = 600

    # Load model and set protocol, create simulation
    m = myokit.load_model('ohara-cipa-v1-2017.mmt')
    p = myokit.pacing.blocktrain(bcl, 0.5, offset=0, level=1.0, limit=0)
    s = myokit.Simulation(m,p)

    # Set cell type
    cell_types = {'Endocardial' : 0, 'Epicardial' : 1, 'Mid-myocardial' : 2}
    cell_type =  'Epicardial'
    s.set_constant('cell.celltype', cell_types[cell_type])

    # Pre-pace simulation to steady state and run 5 cycles of simulation
    #s.pre(400*bcl)

    # Set threshold, 90% of repolarisation to resting membrane potential
    thresh = 0.9 * s.state()[m.get('membrane.V').indice()]

    ss = steady(m,p, pcl = 400, paces = 500)
    print ss
    #s = myokit.Simulation(m,p)
    #d = s.run(paces*bcl, log_interval = 0.1 )


    # Running using function
    #start, duration = ap_duration(m,p, paces = 10,thresh)
    #print start
    #print duration
    pl.show()

if __name__ == "__main__":
    main()
