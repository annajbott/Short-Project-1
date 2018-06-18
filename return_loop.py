#!/usr/bin/env python
import myokit
import matplotlib.pyplot as pl
from manual_APD import ap_duration
import numpy as np


## Return Loop protocol ##

# Plots for protocol, APs and restitution curve
def return_loop(model, number_runs = 25, cell_type = 1, AP_plot = False, protocol_plot = True, restitution_curve = True):

    cell_types = {0:'Endocardial', 1: 'Epicardial', 2: 'Mid-myocardial'}
    if model == 'tentusscher-2006':
        label = 'cell.type'
        name  = 'Ten-Tusscher (2006)'
    elif model == 'grandi-2010':
        label = 'type.epi'
        cell_types = {0:'Endocardial', 1: 'Epicardial'}
        name = 'Grandi (2010)'
    elif model == 'grandi-2010_modified':
        label = 'type.epi'
        cell_types = {0:'Endocardial', 1: 'Epicardial'}
        name = 'Grandi (2010) with late sodium channel'
    elif model == 'ohara-2011':
        label = 'cell.mode'
        name = "O'hara (2011)"
    elif model == 'ohara-cipa-v1-2017.mmt':
        label = 'cell.celltype'
        name = "O'hara- CiPA (2017)"

    m = myokit.load_model('{}.mmt'.format(model))

    # Protocol
    p = myokit.Protocol()


    max_pcl = 600
    min_pcl = 100

    # Make protocol with return_loop included first for no HF (50-90 pattern)

    # Linearly space return loop
    points = np.ndarray.flatten(np.asarray([np.linspace(min_pcl,max_pcl, 50)[:-1], np.linspace(max_pcl,min_pcl, 50)[:-1]]))
    hr_pacing_list = []

    for i in points:

        hr_pacing_list.append(i)

    # Numpy arrays
    hr_pacing_list = np.asarray(hr_pacing_list)

    # Generate integer pacing values from instantaneous HR
    pacing_list = [int(element) for element in hr_pacing_list]
    beats_per_pace = 1

    # Setting initital offset to zero
    offset = 0

    # Empty array to store offsets of new pacing in
    offset_array = []

    # Array to store every PCL in
    pace_array = []

    for number in range(0, number_runs):
        # Create protocol for healthy heart pacing
        for pacing in pacing_list:

            p.schedule(1, start = offset, duration = 0.5, period = pacing, multiplier = beats_per_pace)
            offset += beats_per_pace*pacing
            pace_array.append(pacing)
            offset_array.append(offset)

    ## Create and run simulation for healthy ventricular cell model ##
    ## ------------------------------------------------------------ ##

    s = myokit.Simulation(m, p)
    s.set_constant(label, cell_type)
    s.pre(600*100)
    s.reset()
    d = s.run(offset, log = ['membrane.V','engine.time'])
    start, duration, thresh = ap_duration(d, number_runs*len(pacing_list), repolarisation = 90)

    if AP_plot != False:
        # Plot APs for healthy ventricular cell model
        pl.figure()
        pl.plot(d['engine.time'],d['membrane.V'])
        pl.xlabel('Time (ms)')
        pl.ylabel('Membrane Potential (mV)')

        for i in range(0, 21):
            # Plot arrows and PCL for first sin wave of APs
            pl.arrow(start[i], thresh[i], start[i+1]-start[i], 0, head_width=3, head_length=80,
            length_includes_head=True, color = 'green')
            pl.text(start[i] + 20, -72, str(pace_array[i]))


    if protocol_plot == True:
        # Figure to plot HR (BPM) against time

        pl.figure()
        pl.plot(offset_array,pace_array,'.-')
        pl.xlim(0,60000)
        pl.xlabel('Time (ms)')
        pl.ylabel('PCL (ms)')
        pl.title('Return loop protocol')

    if restitution_curve == True:
        pl.figure()
        if len(duration)< len(pace_array):
            pace_array = pace_array[:-1]
        offset_array = np.asarray(offset_array)
        pcl_start = np.zeros(len(duration))

        pcl_start_up = np.zeros(len(duration))
        pcl_start_down = np.zeros(len(duration))
        duration_up = np.zeros(len(duration))
        duration_down = np.zeros(len(duration))

        length = len(duration)
        j = 0
        q = 0
        for i in range(1,length):
            z = np.nonzero(offset_array < start[i-1])[0][-1]
            #pcl_start[i] = pace_array[z]
            if (z/50)%2 == 0:
                pcl_start_up[j] = pace_array[z]
                duration_up[j] = duration[i]
                j += 1
            else:
                pcl_start_down[q] = pace_array[z]
                duration_down[q] = duration[i]
                q += 1

        #pcl_start = pcl_start[np.nonzero(pcl_start)]
        pcl_start_up = pcl_start_up[np.nonzero(pcl_start_up)]
        pcl_start_down = pcl_start_down[np.nonzero(pcl_start_down)]
        duration_down = duration_down[np.nonzero(duration_down)]
        duration_up = duration_up[np.nonzero(duration_up)]

        pl.plot(pcl_start_up[4*len(pcl_start_up)/5], duration_up[4*len(pcl_start_up)/5], '.')
        pl.plot(pcl_start_down[4*len(pcl_start_down)/5], duration_down[4*len(pcl_start_down)/5], '.')

        pl.xlabel('PCL (ms)')
        pl.ylabel('APD (ms)')
        pl.title('Restitution Curves for Return-loop Protocol')
        pl.legend(['High HR--> Low HR','Low HR--> High HR' ])


    pl.show()
# Main function for testing
def main():
    return_loop(model = 'ohara-2011', cell_type = 0, protocol_plot =  True, restitution_curve = True, AP_plot = False)

if __name__ == "__main__":
    main()
