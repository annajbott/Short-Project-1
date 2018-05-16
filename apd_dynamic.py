#!/usr/bin/env python

import matplotlib.pyplot as pl
import myokit
import numpy as np

def apd_dynamic(d,p, paces, repolarisation = 90):

    ## Get times when paces occur ##
    ## -------------------------- ##

    # Convert membrane potential and time lists to numpy arrays
    V = np.asarray(d['membrane.V'])
    # Depending on model, sometimes environment.time, others engine.time
    time = np.asarray(d['engine.time'])

    pacelog = p.create_log_for_interval(0, p.characteristic_time(), for_drawing = True)
    pace_log = np.asarray(pacelog['pace'])
    pace_log_times =  np.asarray(pacelog.time())

    # Indexes of values where pacing takes place, i.e. 1s not 0s
    pace_index = np.nonzero(pace_log)
    # Uses index numbers to find times in time list, when pacing occurs
    # Removes every second pacing time, those with 0.5 added to the previous value
    pacing_times  = pace_log_times[pace_index][::2]

    ## Get times when PCL changes ##
    ## -------------------------- ##

    pacing_difference = np.zeros(len(pacing_times)-1)
    for i in range(1, len(pacing_times)):
        pacing_difference[i-1] = pacing_times[i] - pacing_times[i-1]

    start_pace_index = np.nonzero(np.ediff1d(pacing_difference, to_begin = 0))

    pace_start = [pacing_times[start_pace_index[i] + 1] for i in range(0,len(start_pace_index))][0]
    pace_start = np.insert(pace_start,0,0)
    # End of simulation at the end
    pace_start = np.insert(pace_start,len(pace_start), time[-1])
    pace_start = np.asarray(pace_start)

    ## Get different thresholds for each PCL region

    pcl_thresh = np.zeros(len(pace_start))
    for i in range(0, len(pace_start) -1):
        pace_start_index = np.nonzero(time >= pace_start[i])[0][0]
        pace_end_index = np.nonzero(time <= pace_start[i+1])[0][-1]
        pace_stable_start_index = pace_start_index + int((pace_end_index - pace_start_index)/3)
        pcl_thresh[i] = min(V[pace_stable_start_index:pace_end_index]) + 5

    #last_pace_stable = np.nonzero(time >= pace_start[-1])[0][0] + int((pace_end_index - pace_start_index)/3)
    #pcl_thresh[-1] = min(V[last_pace_stable:]) + 8


    ## Initialising vales ##
    ## ------------------ ##

    # Blank numpy arrays to contain resting values, max peaks for each AP
    # Times of peaks, duration of APs and start of AP above some threshold
    # Assuming no 1:2 (or 1:3 etc.) ratios, so number of paces >= number of APs
    resting_values = np.zeros(2*paces + 2)
    peak_values = np.zeros(2*paces + 2)
    time_of_peak_values = np.zeros(2*paces + 2)
    onset = np.zeros(2*paces + 2) # Plus 2 in case starts mid AP and ends with incomplete AP
    duration_ap = np.zeros(2*paces + 2)

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


    ## Iterating over time points ##
    ## -------------------------- ##

    for i in range(1, len(V)-1):
        pcl_protocol_number = np.nonzero(time[i] >= pace_start)[0][-1]
        AP_threshold = pcl_thresh[pcl_protocol_number]
        pacing_time_range = int(pace_start[pcl_protocol_number+1]-pace_start[pcl_protocol_number])
        # Only look at APs after things have settled a bit more
        if time[i] < pace_start[pcl_protocol_number] + pacing_time_range/3:
            in_ap = 0
            continue

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

        # If end of pacing protocol for this PCL, end AP and delete its info
        elif pcl_protocol_number != np.nonzero(time[i-1] >= pace_start)[0][-1]:

            # No longer in AP
            in_ap = 0

            # Remove start time
            onset = onset[:-1]

            #Remove resting value
            resting_values = resting_values[:-1]

            # No peak values or times stored, tally not updated.

            # Re-initialise values for next AP
            current_peak_time = -inf
            max_upstroke_velocity = -inf
            current_peak = -inf
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
    thresh = []

    # Iterating over APs rather than time points
    for ap_index in range(0, len(peak_values)):
        # Custom threshold for each AP. Range*(100-repolarisation %) + resting value
        custom_thresh = resting_values[ap_index] + 0.01*(100-repolarisation)*(peak_values[ap_index]-resting_values[ap_index])
        thresh.append(custom_thresh)
        starting_time_index = np.nonzero(time > onset[ap_index])[0][0] - 1
        time_ap = time[starting_time_index]
        pcl_protocol_number = np.nonzero(time_ap >= pace_start)[0][-1]
        AP_threshold = pcl_thresh[pcl_protocol_number]
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
                    #print apd_end_time, apd_start_time
                    break
                prev_t = time[t]
                prev_v = V[t]

        # Re-initialise for next AP.
        apd_starting_index = None
        apd_end_time = None
        apd_start_time = None
        # Loop to next AP

    duration_ap = duration_ap[np.nonzero(duration_ap)]
    return[onset_apd, duration_ap, thresh]
