#!/usr/bin/env python
import myokit
import matplotlib.pyplot as pl
from manual_APD import ap_duration
import numpy as np
from HF_model import *

## Return Loop protocol ##

# Plots for protocol, APs and restitution curve
def return_loop(model, HF_model, number_runs = 25, cell_type = 1, AP_plot = False, protocol_plot = True, restitution_curve = True, APD_time_plot = False):

    cell_types = {0:'Endocardial', 1: 'Epicardial', 2: 'Mid-myocardial'}
    if model == 'tentusscher-2006':
        label = 'cell.type'
        name  = 'Ten-Tusscher (2006)'
        HF_label = 'TT'
    elif model == 'grandi-2010':
        label = 'type.epi'
        cell_types = {0:'Endocardial', 1: 'Epicardial'}
        name = 'Grandi (2010)'
        HF_label = 'GPB'
    elif model == 'grandi-2010_modified':
        label = 'type.epi'
        cell_types = {0:'Endocardial', 1: 'Epicardial'}
        name = 'Grandi (2010) with late sodium channel'
        HF_label = 'GPB'
    elif model == 'ohara-2011':
        label = 'cell.mode'
        name = "O'hara (2011)"
        HF_label = 'Ord'
    elif model == 'ohara-cipa-v1-2017.mmt':
        label = 'cell.celltype'
        name = "O'hara- CiPA (2017)"
        HF_label = 'Ordcipa'

    m = myokit.load_model('{}.mmt'.format(model))
    models_HF = {'Ord_HF_Gomez' : Ord_HF_Gomez, 'Ord_HF_Elshrif' : Ord_HF_Elshrif, 'GPB_HF_Gomez' : GPB_HF_Gomez, 'GPB_HF_Gomez': GPB_HF_Moreno, 'TT_HF_Lu' : TT_HF_Lu}
    if HF_model != None:
        m_str = '{}_HF_{}'.format(HF_label, HF_model)
        m_func = models_HF[m_str]
        hf_m = m_func(cell_type)


    # Protocol
    p = myokit.Protocol()
    hf_p = myokit.Protocol()


    max_pcl = 600
    min_pcl = 150

    # Make protocol with return_loop included first for no HF (50-90 pattern)

    # Linearly space return loop
    points = np.ndarray.flatten(np.asarray([np.linspace(600,100, 50)[:-1], np.linspace(100,600, 50)[:-1]]))
    hr_pacing_list = []
    hf_hr_pacing_list = []

    for i in points:

        hr_pacing_list.append(i)
        hf_hr_pacing_list.append(i)

    # Numpy arrays
    hr_pacing_list = np.asarray(hr_pacing_list)
    hf_hr_pacing_list = np.asarray(hf_hr_pacing_list)

    # Generate integer pacing values from instantaneous HR
    pacing_list = [int(element) for element in hr_pacing_list]
    hf_pacing_list = [int(element)  for element in hf_hr_pacing_list]
    beats_per_pace = 1

    # Setting initital offset to zero
    offset = 0
    hf_offset = 0


    # Empty array to store offsets of new pacing in
    offset_array = []
    hf_offset_array = []

    # Array to store every PCL in
    pace_array = []
    hf_pace_array = []

    for number in range(0, number_runs):
        # Create protocol for healthy heart pacing
        for pacing in pacing_list:

            p.schedule(1, start = offset, duration = 0.5, period = pacing, multiplier = beats_per_pace)
            offset += beats_per_pace*pacing
            pace_array.append(pacing)
            offset_array.append(offset)
        # Create protocol for HF at rest
        for pacing in hf_pacing_list:
            hf_p.schedule(1, start = hf_offset, duration = 0.5, period = pacing, multiplier = beats_per_pace)
            hf_offset += beats_per_pace*pacing
            hf_pace_array.append(pacing)
            hf_offset_array.append(hf_offset)


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

    ## Create and run simulation for HF at rest ##
    ## ---------------------------------------- ##

    hf_s = myokit.Simulation(hf_m, hf_p)
    hf_s.set_constant(label, cell_type)
    hf_s.pre(np.sum(hf_pacing_list)*10)
    hf_s.reset()
    hf_d = hf_s.run(hf_offset, log = ['membrane.V','engine.time'])
    # Calculate APDs
    hf_start, hf_duration, hf_thresh = ap_duration(hf_d, number_runs*len(pacing_list), repolarisation = 90)


    if AP_plot != False:
        # Plot APs for HF ventricular cell model at rest
        pl.figure()
        pl.plot(hf_d['engine.time'],hf_d['membrane.V'])
        pl.title('HF')

        for i in range(0, 21):
            # Plot arrows and PCL for first sin wave of APs
            pl.arrow(hf_start[i], hf_thresh[i], hf_start[i+1]-hf_start[i], 0, head_width=3, head_length=80,
            length_includes_head=True, color = 'green')
            pl.text(hf_start[i] + 20, -72, str(hf_pace_array[i]))

    if protocol_plot == True:
        # Figure to plot HR (BPM) against time
        hr_bpm = [60*1000.0/element for element in pace_array]
        hf_hr_bpm = [60*1000.0/element for element in hf_pace_array]

        pl.figure()
        pl.plot(offset_array,hr_bpm,'.-')
        pl.xlim(0,50000)
        pl.xlabel('Time (ms)')
        pl.ylabel('Instantaneous HR (BPM)')
        pl.title('Return loop protocol')

    if restitution_curve == True:
        pl.figure()
        if len(duration)< len(pace_array):
            pace_array = pace_array[:-1]
        offset_array = np.asarray(offset_array)
        pcl_start = np.zeros(len(duration))
        hf_pcl_start = np.zeros(len(hf_duration))

        pcl_start_up = np.zeros(len(duration))
        pcl_start_down = np.zeros(len(duration))
        duration_up = np.zeros(len(duration))
        duration_down = np.zeros(len(duration))



        length = max(len(duration), len(hf_duration))
        j = 0
        q = 0
        for i in range(0,length):
            z = np.nonzero(offset_array < start[i])[0][-1]
            pcl_start[i] = pace_array[z]
            if (z/50)%2 == 0:
                pcl_start_up[j] = pace_array[z]
                duration_up[j] = duration[i]
                j += 1
            else:
                pcl_start_down[q] = pace_array[z]
                duration_down[q] = duration[i]
                q += 1

        #    hf_z = np.nonzero(hf_offset_array < hf_start[i])[0][-1]
        #    hf_pcl_start[i] = hf_pace_array[hf_z]
        pcl_start = pcl_start[np.nonzero(pcl_start)]
        pcl_start_up = pcl_start_up[np.nonzero(pcl_start_up)]
        pcl_start_down = pcl_start_down[np.nonzero(pcl_start_down)]
        duration_down = duration_down[np.nonzero(duration_down)]
        duration_up = duration_up[np.nonzero(duration_up)]

        #hf_pcl_start = hf_pcl_start[np.nonzero(hf_pcl_start)]
        #pl.plot(pcl_start[len(pcl_start)/2 -1:-1], duration[len(pcl_start)/2:], '.')
        pl.plot(pcl_start_up[-100:-1], duration_up[-99:], '.')
        pl.plot(pcl_start_down[-100:-1], duration_down[-99:], '.')

        #pl.plot(hf_pace_array[len(pace_array)/2 -1: -1], hf_duration[len(pace_array)/2 :], '.')
        pl.xlabel('PCL (ms)')
        pl.ylabel('APD (ms)')
        pl.title('Restitution Curves for Return-loop Protocol')
        pl.legend(['Low HR--> High HR', 'High HR--> Low HR'])

    ## APD vs time plot ##
    ## ---------------- ##

    if APD_time_plot == True:
        pl.figure()
        pl.plot(offset_array, duration,'.-')
        pl.plot(hf_offset_array, hf_duration,'.-')
        pl.xlabel('Time (ms)')
        pl.ylabel('APD (ms)')
        pl.title('APD varying over time with the return_loop protocol')
    # Show plots
    pl.show()


# Main function for testing
def main():
    return_loop(model = 'ohara-2011', HF_model = 'Gomez', cell_type = 0, protocol_plot =  False, APD_time_plot = False, restitution_curve = True, AP_plot = False)

if __name__ == "__main__":
    main()