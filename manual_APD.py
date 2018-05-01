#!/usr/bin/env python

import matplotlib.pyplot as pl
import myokit
import numpy as np

def ap_duration(d, paces, repolarisation = 90 , time_interval = None):

    # Convert membrane potential and time lists to numpy arrays
    V = np.asarray(d['membrane.V'])
    # Depending on model, sometimes environment.time, others engine.time
    time = np.asarray(d['engine.time'])

    # Blank numpy arrays to contain resting values, max peaks for each AP
    # Times of peaks, duration of APs and start of AP above some threshold
    # Assuming no 1:2 (or 1:3 etc.) ratios, so number of paces >= number of APs
    resting_values = np.zeros(paces + 2)
    peak_values = np.zeros(paces + 2)
    time_of_peak_values = np.zeros(paces + 2)
    onset = np.zeros(paces + 2) # Plus 2 in case starts mid AP and ends with incomplete AP
    duration_ap = np.zeros(paces + 2)

    # Segmenting each period using the minimum voltage + 5mV as threshold
    AP_threshold = min(V) + 5
    # Boolean whether in AP or outside
    in_ap = 0
    # Counter for number of APs
    responses = 0
    max_upstroke_velocity = 0.0 # Unsure if needed for this purpose
    inf = float("inf")

    current_minimum_velocity = inf
    current_resting_value = inf
    current_peak = -inf
    resting_potential_gradient_thresh = 0.01 # Gradient taken from Chaste code

    # Boolean for if section of AP where gradient < 0.01 has been found
    flat_bit_found = 0
    apd_end_time = None
    apd_start_time = None

    # Iterating over time points
    for i in range(1,len(V)):
        # Taken from Chaste code: using gradients to determine upstrokes and rest

        # Gradient calculated using difference in potential over time difference from previous time point
        voltage_grad = (V[i]-V[i-1])/float((time[i]-time[i-1]))
        if voltage_grad >= max_upstroke_velocity:
            max_upstroke_velocity = voltage_grad
            current_time_upstroke = time[i]

        # Looking for rest section of AP, either flat-ish gradient or minimum voltage
        if abs(voltage_grad) <= current_minimum_velocity and abs(voltage_grad) <= resting_potential_gradient_thresh and not in_ap:
            current_minimum_velocity = abs(voltage_grad)
            current_resting_value = V[i-1]
            flat_bit_found  = 1

        # Keep storing mimimum voltage if flat bit has not been found for period
        elif V[i-1] < current_resting_value and flat_bit_found == 0 and not in_ap:
            current_resting_value = V[i-1]

        ### Crossing threshold. Start of AP ###
        if V[i] > AP_threshold and V[i-1] <= AP_threshold and not in_ap:

            # Now within AP
            in_ap = 1

            # Set resting value for this AP. Used as minimum
            resting_values[responses] = current_resting_value

            # Re-initialise these values for next AP
            current_minimum_velocity = inf
            current_resting_value = inf
            flat_bit_found = 0

            # Linear interpolation to store onset time.
            # y = y0 + (x -x0)*(y1-y0)/(x1 -x0). y: time, x : voltage
            onset[responses] = time[i-1] + (AP_threshold - V[i-1])*(time[i]-time[i-1])/(V[i]-V[i-1])

            switching_phase = True

        # In AP, record height of the peak and corresponding time. Useful check if printed
        elif in_ap and V[i] > current_peak:
            current_peak = V[i]
            current_peak_time = time[i]

        # Cross threshold again. End of AP. Update variables
        elif in_ap and V[i] < AP_threshold and V[i-1] >= AP_threshold:

            # No longer in AP
            in_ap = 0

            # Record peak values
            peak_values[responses] = current_peak
            time_of_peak_values[responses] = current_peak_time

            # Updates tally of number of APs
            responses += 1

            # Re-initialise values for next AP
            current_peak_time = -inf
            max_upstroke_velocity = -inf
            current_peak = AP_threshold
            current_time_of_upstroke_velocity = 0.0

        # Loop round to next time point

    # Keep only non-zero terms. Matrices containing start and duration of APs at set threshold
    onset = onset[np.nonzero(onset)]
    peak_values = peak_values[np.nonzero(peak_values)]
    time_of_peak_values = time_of_peak_values[np.nonzero(time_of_peak_values)]

    ### Varying thresholds calculated for each AP ###
    ### ----------------------------------------- ###
    # Using min and max of each AP, calculating specific percent of AP duration

    # Clone onset matrix, to fill with APD (90/50 etc) starts
    onset_apd = onset

    # Iterating over APs rather than time points
    for ap_index in range(0, len(peak_values)):
        # Custom threshold for each AP. Range*(100-repolarisation %) + resting value
        custom_thresh = resting_values[ap_index] + 0.01*(100-repolarisation)*(peak_values[ap_index]-resting_values[ap_index])

        starting_time_index = np.nonzero(time > onset[ap_index])[0][0] - 1

        # The APD 90/50.. threshold is greater than earlier threshold
        # Look forwards in time
        if custom_thresh >= AP_threshold:
            prev_v = V[starting_time_index]
            prev_t = time[starting_time_index]
            for t in range(starting_time_index,len(time)):
                if V[t] > custom_thresh and V[t-1] < custom_thresh:
                    # Linear interpolation
                    apd_start_time = prev_t + ((custom_thresh - prev_v) / float(V[t] - prev_v)) * (time[t] - prev_t)
                    onset_apd[ap_index] = apd_start_time
                    apd_starting_index = t
                    break
                prev_v = V[t]
                prev_t = time[t]

        # Look backwards to find custom threshold, i.e. threshold is lower than previous threshold used
        else:
            prev_v = V[starting_time_index + 1]
            prev_t = time[starting_time_index + 1]
            for t in range(starting_time_index, 0, -1):
                if V[t] < custom_thresh:
                    # Linear interpolation
                    apd_start_time = prev_t + ((custom_thresh - prev_v) / float((V[t] - prev_v))) * (time[t] - prev_t)
                    onset_apd[ap_index] = apd_start_time
                    apd_starting_index = (t + 1);
                    break;
                prev_v = V[t]
                prev_t = time[t]

        if apd_end_time != None and apd_start_time != None and apd_start_time < apd_end_time:
            # Skip to next AP if start of this AP is reached before last AP is finished
            continue


        # If apd_start_index has a value then look forwards in time for repolarisation
        if apd_starting_index != None:
            prev_t = time[apd_starting_index - 1]
            prev_v = V[apd_starting_index -1]
            for t in range(apd_starting_index, len(time)):
                if V[t] < custom_thresh:
                    # Linear interpolation
                    apd_end_time = time[t - 1] + ((custom_thresh - prev_v) / (V[t] - prev_v)) * (time[t] - time[t - 1])
                    duration_ap[ap_index] = apd_end_time - apd_start_time
                    break
                prev_t = time[t]
                prev_v = V[t]

        # Re-initialise for next AP.
        apd_starting_index = None
        apd_end_time = None
        apd_start_time = None
        # Loop to next AP

    duration_ap = duration_ap[np.nonzero(duration_ap)]
    return[onset_apd, duration_ap]

def steady(m,p,pcl,paces):
    # Create simulation using user defined model and protocol
    s = myokit.Simulation(m,p)

    # Progress report
    print "In steady function"

    # Run for pcl*paces, specifying to use time values every 0.1 ms interval for solver
    d  = s.run(paces*pcl, log = ['engine.time', 'membrane.V'], log_interval = 0.1)

    # Lengthy step due to log_interval component, so progress report when this is finished
    print "Simulations finished running"

    # Use function to calculate start and AP duration, default APD 90
    start, duration = ap_duration(d, paces)

    # Find long-short patterns. Linearly repeating
    for i in range(1, 10):
        if np.round(duration[-1],0) != np.round(duration[-1 - i],0):
            continue
        else:
            seq = i
            print seq
            break

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
    perc = np.zeros(len(V)-seq)
    dist_array = np.zeros(len(V)-seq)

    # Length V = number of rows --> number of APs (periods)
    for i in range(seq,len(V)):
        # Euclidean dist
        #dist = np.sqrt(np.sum((V[i]-V[i-seq])**2))

        # L-infinite
        dist = max(abs(V[i]-V[i-seq]))


        #perc[i-seq] = (dist/abs(np.mean(V[i-seq])))*100
        dist_array[i- seq] = dist
    # Moving average for 10 terms calculated for percentage change
    #moving_av = np.convolve(perc, np.ones((10,))/10, mode='valid')
    moving_av_dist = np.convolve(dist_array, np.ones((10,))/10, mode='valid')
    pl.figure()
    pl.plot(d['engine.time'],d['membrane.V'])

    # Plot results
    pl.figure()
    pl.plot(range(seq,len(V)), dist_array)
    pl.xlabel('Period Number')
    pl.ylabel("L-infinite norm of membrane potential and previous period's membrane potential")

    if np.size(np.nonzero(moving_av_dist < 0.5)) < 1:
        ss =  "Run again with more paces"
    else:
        ss = paces_per_period*np.nonzero(moving_av_dist < 0.4)[0][0]

    return(ss)

def main():
    ## Initialisation common to steady and APD calc ##
    # Set pacing
    bcl = 150

    # Load model and set protocol, create simulation
    m = myokit.load_model('ohara-cipa-v1-2017.mmt')
    p = myokit.pacing.blocktrain(bcl, 0.5, offset=0, level=1.0, limit=0)
    s = myokit.Simulation(m,p)

    # Set cell type
    cell_types = {'Endocardial' : 0, 'Epicardial' : 1, 'Mid-myocardial' : 2}
    cell_type =  'Epicardial'
    s.set_constant('cell.celltype', cell_types[cell_type])


    ### Testing APD calc section ###
    ### ------------------------ ###
    # Pre-pace simulation and run 30 cycles of simulation
    #'''
    s.pre(200*bcl)
    # Running using function
    paces = 10
    d = s.run(paces*bcl)
    start, duration = ap_duration(d, paces)
    print start
    print duration
    pl.figure()
    pl.plot(d['engine.time'], d['membrane.V'])

    '''
    ### Testing steady state calc ###
    ### ------------------------- ###
    ss = steady(m,p, bcl, paces = 300)
    print ss
    '''
    pl.show()
if __name__ == "__main__":
    main()
