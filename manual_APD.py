#!/usr/bin/env python

import matplotlib.pyplot as pl
import myokit
import numpy as np

def ap_duration(d, paces, repolarisation = 90 , time_interval = None):

    # Convert membrane potential and time arrays to numpy
    V = np.asarray(d['membrane.V'])
    time = np.asarray(d['environment.time'])

    # Numpy arrays containing start of APs and duration of APs
    # Assuming no 1:2 (or 1:3 etc) ratios, so number of paces >= number of APs
    start_ap = np.zeros(paces + 2) # Plus 2 in case starts mid AP and ends with incomplete AP
    start_indexes = np.zeros(paces + 2)
    duration_ap = np.zeros(paces + 2)
    end_ap = np.zeros(paces + 2)
    end_indexes = np.zeros(paces + 2)

    # Segmenting each period using -75mV as AP threshold.
    AP_threshold = -80
    in_ap = 0
    responses = 0
    for i in range(0,len(V)):
        # Start of AP (so not already in AP), greater than threshold
        if V[i] >= AP_threshold and not in_ap:
            in_ap = 1
            start_ap[responses] = time[i]
            start_indexes[responses] = i
        # End of AP, record end time and subtract start for the duration
        elif in_ap and V[i] < AP_threshold:
            in_ap = 0
            end_ap[responses] = time[i]
            end_indexes[responses] = i
            duration_ap[responses] = end_ap[responses] - start_ap[responses]
            responses += 1

    # If simulation was started in the middle of an AP, discard first entry
    if V[0] > AP_threshold:
        start_ap = start_ap[1:]
        duration_ap = duration_ap[1:]
    # If final AP not finished, discard its start time
    if V[-1] > AP_threshold:
        start_ap = start_ap[0:-1]

    # Keep only non-zero terms. Matrices containing start and duration of APs at set threshold
    start_ap = start_ap[np.nonzero(start_ap)]
    end_ap = end_ap[np.nonzero(end_ap)]
    duration_ap = duration_ap[np.nonzero(duration_ap)]

    ### Varying thresholds calculated for each AP ###
    ### ----------------------------------------- ###
    # 2nd part if want the thresholds to vary for the size of the AP

    # Clone matrices from earlier section
    duration2_ap = duration_ap
    start2_ap = start_ap

    for i in range(0,len(duration2_ap)):
        #print i
        if i == 0:
            # For first action potential
            ap_period_start = 0
            ap_period_end = int(start_indexes[i+1])

        elif i == len(start_ap) - 1:
            # For middle APs
            ap_period_start = int(end_indexes[i-1])
            ap_period_end = None

        else:
            # For final action potential
            ap_period_start = int(end_indexes[i-1]) + 1
            ap_period_end  = int(start_indexes[i+1]) -1

        ap_period = V[ap_period_start:ap_period_end]

        # Calculate range and threhsold for each individual AP
        custom_range = (max(ap_period) - min(ap_period))
        custom_thresh = min(ap_period) + custom_range*(100-repolarisation)*0.01

        # Cut off the front end and back end of each interval, ensure not exceeding threshold already.
        while ap_period[0] >= custom_thresh:
            ap_period = ap_period[1:]

        while ap_period[-1] >= custom_thresh:
            ap_period = ap_period[:-1]

        # Find the first and last instance where threshold is exceeded.
        start_index = np.nonzero(ap_period> custom_thresh)[0][0] + ap_period_start
        end_index = np.nonzero(ap_period> custom_thresh)[0][-1] + ap_period_start
        duration2_ap[i] = time[end_index] - time[start_index]
        start2_ap[i] = time[start_index]


    return[start2_ap, duration2_ap]

def steady(m,p,pcl,paces):
    # Create simulation using user defined model and protocol
    s = myokit.Simulation(m,p)

    # Progress report
    print "In steady function"

    # Run for pcl*paces, specifying to use time values every 0.1 ms interval for solver
    d  = s.run(paces*pcl, log_interval = 0.1)

    # Lengthy step due to log_interval component, so progress report when this is finished
    print "Simulations finished running"

    # Use function to calculate start and AP duration, threshold -75mV
    start, duration = ap_duration(d, paces, threshold = -75)


    # If there is greater than a 1:1 ratio
    paces_per_period = np.round(float(start[-1] - start[-2]),0)/ float(pcl)

    # How many time points in each period (AP) (time interval for solver 0.1)
    time_points_per_period = pcl*paces_per_period/0.1

    # Convert membrane potential list into numpy array
    V = np.asarray(d['membrane.V'])
    # Long-short pattern
    if np.round(duration[-1],0) != np.round(duration[-2],0) and np.round(duration[-1],0) == np.round(duration[-3],0):
        long_short = 1
    else:
        long_short = 0

    # If multiple paces per period, V may not be divisible by no. time_points_per_period
    # Slice the end of the array off so it is divisble and can reshape array
    remainder = int(len(V)%time_points_per_period)

    if remainder != 0:
        V = V[0:-remainder]

    # Reshape array so that each period is on a seperate row
    # Each column contains the same point in an AP, but for different beats
    V = np.reshape(V, (-1, int(time_points_per_period)))

    # Euclidean distance and percentage change between adjacent periods calculated
    perc = np.zeros(len(V)-1 - long_short)
    dist_array = np.zeros(len(V)-1 - long_short)

    # Length V = number of rows --> number of APs (periods)
    for i in range(1 + long_short,len(V)):
        dist = np.sqrt(np.sum((V[i]-V[i-1 - long_short])**2))
        #perc[i-1 - long_short] = (dist/abs(np.mean(V[i-1 - long_short])))*100
        dist_array[i-1 - long_short] = dist
    # Moving average for 10 terms calculated for percentage change
    #moving_av = np.convolve(perc, np.ones((10,))/10, mode='valid')
    moving_av_dist = np.convolve(dist_array, np.ones((10,))/10, mode='valid')
    pl.figure()
    pl.plot(d['environment.time'],d['membrane.V'])

    # Plot results
    pl.figure()
    pl.plot(range(1 + long_short,len(V)), dist_array)
    pl.xlabel('Period Number')
    pl.ylabel("Euclidean distance of membrane potential from previous period's membrane potential")
    '''
    pl.figure()
    pl.plot(range(1 + long_short,len(V)-9), moving_av)
    pl.xlabel('Period Number')
    pl.ylabel('Moving average')

    if np.size(np.nonzero(moving_av < 1)) < 1:
        ss =  "Run again with more paces"
    else:
        ss = paces_per_period*np.nonzero(moving_av < 1)[0][0]
    '''
    if np.size(np.nonzero(moving_av_dist < 0.8)) < 1:
        ss =  "Run again with more paces"
    else:
        ss = paces_per_period*np.nonzero(moving_av_dist < 0.8)[0][0]
    return(ss)

def main():
    #Set pacing
    bcl = 80

    # Load model and set protocol, create simulation
    m = myokit.load_model('ohara-cipa-v1-2017.mmt')
    p = myokit.pacing.blocktrain(bcl, 0.5, offset=0, level=1.0, limit=0)
    s = myokit.Simulation(m,p)

    # Set cell type
    cell_types = {'Endocardial' : 0, 'Epicardial' : 1, 'Mid-myocardial' : 2}
    cell_type =  'Epicardial'
    s.set_constant('cell.celltype', cell_types[cell_type])

    # Pre-pace simulation to steady state and run 5 cycles of simulation
    s.pre(400*bcl)

    # Set threshold, 90% of repolarisation to resting membrane potential
    # thresh = 0.9 * s.state()[m.get('membrane.V').indice()]

    paces = 30
    #ss = steady(m,p, bcl, paces = 1000)
    #print ss
    #s = myokit.Simulation(m,p)
    d = s.run(paces*bcl)#, log_interval = 0.1 )


    # Running using function
    start, duration = ap_duration(d, paces = 30)
    print start
    print duration
    pl.figure()
    pl.plot(d['environment.time'], d['membrane.V'])
    pl.show()

if __name__ == "__main__":
    main()
