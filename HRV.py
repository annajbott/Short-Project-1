#!/usr/bin/env python
import myokit
import matplotlib.pyplot as pl
from manual_APD import ap_duration
import numpy as np
from HF_model import *


## Heart Rate Variability (HRV)- Healthy, HF (rest), HF (submaximal exercise) ##
## -------------------------------------------------------------------------- ##

# Plots for protocol, APs and restitution curve
def HRV(model, HF_model, number_runs = 50, cell_type = 1, noise = False, AP_plot = False, protocol_plot = True, restitution_curve = True, APD_time_plot = False):

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
    hf_6minwalk_p = myokit.Protocol()


    # 0.15-0.4Hz breathing rate (9-24 breaths a minute / 2500ms - 6666ms)
    # Make protocol with HRV included first for no HF (50-90 pattern)

    # Using sin wave to simulate step protocol
    points = np.linspace(0,2*np.pi, num  = 30, endpoint = False)
    hr_pacing_list = []
    hf_hr_pacing_list = []
    hf_6minwalk_pacing_list = []
    average = 70
    hf_average = 90
    # Weighted average from Correa, Silva Alves, Bianchim (2013)
    hf_6minwalk_average = 114
    for i in points:
        # Using a sin wave centered around the average (18.5,10, 5 = amplitude)
        hr_pacing_list.append(average +18.5*np.sin(i-np.pi/2.0))
        hf_hr_pacing_list.append(hf_average +10*np.sin(i-np.pi/2.0))
        hf_6minwalk_pacing_list.append(hf_6minwalk_average +8*np.sin(i-np.pi/2.0))
    '''
    if noise == True:
        # Introducing fluctuations around pacing values
        hr_pacing_list = [element + np.random.normal(0, 1) for element in hr_pacing_list]
        hf_hr_pacing_list = [element + np.random.normal(0, 1)  for element in hf_hr_pacing_list]
        hf_6minwalk_pacing_list = [element + np.random.normal(0, 1)  for element in hf_6minwalk_pacing_list]
    '''
    # Numpy arrays
    hr_pacing_list = np.asarray(hr_pacing_list)
    hf_hr_pacing_list = np.asarray(hf_hr_pacing_list)
    hf_6minwalk_pacing_list = np.asarray(hf_6minwalk_pacing_list)

    # Generate integer pacing values from instantaneous HR
    pacing_list = [np.round(60*1000.0/element,0) for element in hr_pacing_list]
    hf_pacing_list = [np.round(60*1000.0/element,0)  for element in hf_hr_pacing_list]
    exercise_pacing_list = [np.round(60*1000.0/element,0)  for element in hf_6minwalk_pacing_list]
    beats_per_pace = 1

    # Setting initital offset to zero
    offset = 0
    hf_offset = 0
    hf_6minwalk_offset = 0

    # Empty array to store offsets of new pacing in
    offset_array = []
    hf_offset_array = []
    hf_6minwalk_offset_array =[]

    # Array to store every PCL in
    pace_array = []
    hf_pace_array = []
    hf_6minwalk_pace_array = []

    for number in range(0, number_runs):
        # Create protocol for healthy heart pacing
        for pacing in pacing_list:
            if noise == True:
                pacing = np.round(np.random.normal(0, 5) + pacing,0)
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
        # Create protocol for HF during final minutes of 6-min walk test
        for pacing in exercise_pacing_list:
            hf_6minwalk_p.schedule(1, start = hf_6minwalk_offset, duration = 0.5, period = pacing, multiplier = beats_per_pace)
            hf_6minwalk_offset += beats_per_pace*pacing
            hf_6minwalk_pace_array.append(pacing)
            hf_6minwalk_offset_array.append(hf_6minwalk_offset)


    ## Create and run simulation for healthy ventricular cell model ##
    ## ------------------------------------------------------------ ##

    s = myokit.Simulation(m, p)
    s.set_constant(label, cell_type)
    s.pre(np.sum(pacing_list)*10)
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


    ## Create and run simulation for HF during final minutes of 6 min walk test ##
    ## ------------------------------------------------------------------------ ##

    hf_6minwalk_s = myokit.Simulation(hf_m, hf_6minwalk_p)
    hf_6minwalk_s.set_constant(label, cell_type)
    hf_6minwalk_s.pre(np.sum(exercise_pacing_list)*10)
    hf_6minwalk_s.reset()
    hf_6minwalk_d = hf_6minwalk_s.run(hf_6minwalk_offset, log = ['membrane.V','engine.time'])
    hf_6minwalk_start, hf_6minwalk_duration, hf_6minwalk_thresh = ap_duration(hf_6minwalk_d, number_runs*len(pacing_list), repolarisation = 90)

    if AP_plot != False:
        # Plot APs for HF ventricular cell model 6 min walk
        pl.figure()
        pl.plot(hf_6minwalk_d['engine.time'],hf_6minwalk_d['membrane.V'])
        pl.title('HF 6-min walk test')

        for i in range(0, 21):
            # Plot arrows and PCL for first sin wave of APs
            pl.arrow(hf_6minwalk_start[i], hf_6minwalk_thresh[i], hf_6minwalk_start[i+1]-hf_6minwalk_start[i], 0, head_width=3, head_length=80,
            length_includes_head=True, color = 'green')
            pl.text(hf_6minwalk_start[i] + 20, -72, str(hf_6minwalk_pace_array[i]))
            if  i == 20:
                pl.arrow(hf_6minwalk_start[i], hf_6minwalk_thresh[i] - 15, hf_6minwalk_duration[i], 0, head_width=3, head_length=80,
                length_includes_head=True, color = 'red')
                pl.text(hf_6minwalk_start[i] + 20, -80, str(hf_6minwalk_duration[i]))


    if protocol_plot == True:
        # Figure to plot HR (BPM) against time
        hr_bpm = [60*1000.0/element for element in pace_array]
        hf_hr_bpm = [60*1000.0/element for element in hf_pace_array]
        hf_6minwalk_bpm = [60*1000.0/element for element in hf_6minwalk_pace_array]

        pl.figure()
        pl.plot(offset_array,hr_bpm,'.-')
        pl.plot(hf_offset_array, hf_hr_bpm,'.-')
        pl.plot(hf_6minwalk_offset_array, hf_6minwalk_bpm,'.-')
        pl.xlim(0,50000)
        pl.xlabel('Time (ms)')
        pl.ylabel('Instantaneous HR (BPM)')
        pl.title('Dynamic Pacing Protocol to mimic HRV')
        pl.legend(['Healthy HRV', 'HF HRV at rest', 'HF HRV 6-min walk test'])

    if restitution_curve == True:
        pl.figure()
        if len(duration)< len(pace_array):
            pace_array = pace_array[:-1]
        pl.plot(pace_array[len(pace_array)/2 -1:-1], duration[len(pace_array)/2:], '.')
        pl.plot(hf_pace_array[len(pace_array)/2 -1: -1], hf_duration[len(pace_array)/2 :], '.')
        pl.plot(hf_6minwalk_pace_array[len(pace_array)/2 -1:-1], hf_6minwalk_duration[len(pace_array)/2:], '.')
        pl.xlabel('PCL (ms)')
        pl.ylabel('APD (ms)')
        pl.title('Restitution Curves for Dynamic Protocol based on HRV')
        pl.legend(['Healthy', 'HF at rest', 'HF 6-min walk test'])

    ## APD vs time plot ##
    ## ---------------- ##

    if APD_time_plot == True:
        pl.figure()
        pl.plot(offset_array, duration,'.-')
        pl.plot(hf_offset_array, hf_duration,'.-')
        pl.plot(hf_6minwalk_offset_array, hf_6minwalk_duration,'.-')
        pl.xlabel('Time (ms)')
        pl.ylabel('APD (ms)')
        pl.title('APD varying over time with the HRV protocol')
    # Show plots
    pl.show()


# Main function for testing
def main():
    HRV(model = 'ohara-2011', HF_model = 'Gomez', noise = False, cell_type = 0, APD_time_plot = True)

if __name__ == "__main__":
    main()
