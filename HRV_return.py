#!/usr/bin/env python
import myokit
import matplotlib.pyplot as pl
from manual_APD import ap_duration
import numpy as np
from HF_model import *


## Heart Rate Variability (HRV)- Healthy, HF (rest), HF (submaximal exercise) ##
## -------------------------------------------------------------------------- ##

# Plots for protocol, APs and restitution curve
def HRV_return(model, HF_model = None, HF_protocol = None, number_runs = 50, cell_type = 0, AP_plot = False, protocol_plot = True, max_PCL = None, min_PCL = None, APD_time_plot = False):

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
    elif model == 'ohara-cipa-v1-2017':
        label = 'cell.celltype'
        name = "O'hara- CiPA (2017)"
        HF_label = 'Ordcipa'

    m = myokit.load_model('{}.mmt'.format(model))
    models_HF = {'Ord_HF_Gomez' : Ord_HF_Gomez, 'Ord_HF_Elshrif' : Ord_HF_Elshrif, 'GPB_HF_Gomez' : GPB_HF_Gomez, 'GPB_HF_Gomez': GPB_HF_Moreno, 'TT_HF_Lu' : TT_HF_Lu,'Ordcipa_HF_Elshrif' : Ordcipa_HF_Elshrif, 'Ordcipa_HF_Gomez': Ordcipa_HF_Gomez}
    if HF_model != None:
        m_str = '{}_HF_{}'.format(HF_label, HF_model)
        m_func = models_HF[m_str]
        m = m_func(cell_type)
        HF_model_label = 'HF ' + HF_model


    # Protocol
    p = myokit.Protocol()

    # Using sin wave to simulate step protocol
    points = np.linspace(0,2*np.pi, num  = 30, endpoint = False)
    hr_pacing_list = []
    if HF_protocol == None:
        average = 75
        amplitude = 19
        colour = 'blue'
        protocol_label = ' Healthy HRV Protocol'
    elif HF_protocol == 'HF':
        average = 90
        amplitude = 12
        protocol_label = ' HF Rest HRV Protocol'

    elif HF_protocol == 'HF_walk':
        # Weighted average from Correa, Silva Alves, Bianchim (2013)
        average = 120
        amplitude = 8
        protocol_label = ' HF 6-min Walk HRV Protocol'

    # If min and max pcl specified, override protocol options
    if max_PCL != None and min_PCL != None:
        max_hr = 60000.0/max_PCL
        min_hr = 60000.0/min_PCL
        average = np.round((max_hr + min_hr)/2.0,0)
        amplitude = np.round((max_hr - min_hr)/2.0,0)
    for i in points:
        # Using a sin wave centered around the average (18.5,10, 5 = amplitude)
        hr_pacing_list.append(average +amplitude*np.sin(i-np.pi/2.0))

    # Numpy arrays
    hr_pacing_list = np.asarray(hr_pacing_list)

    # Generate integer pacing values from instantaneous HR
    pacing_list = [np.round(60*1000.0/element,0) for element in hr_pacing_list]
    beats_per_pace = 1

    # Setting initital offset to zero
    offset = 0

    # Empty array to store offsets of new pacing in
    offset_array = []

    # Array to store every PCL in
    pace_array = []

    for number in range(0, number_runs):
        for pacing in pacing_list:
            # Create protocol schedule
            p.schedule(1, start = offset, duration = 0.5, period = pacing, multiplier = beats_per_pace)
            offset += beats_per_pace*pacing
            pace_array.append(pacing)
            offset_array.append(offset)

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

    offset_array_sec = [element/1000.0 for element in offset_array]

    if protocol_plot == True:
        # Figure to plot HR (BPM) against time
        hr_bpm = [60*1000.0/element for element in pace_array]

        pl.figure()
        pl.plot(offset_array_sec,hr_bpm,'.-')
        pl.xlim(0,50)
        pl.xlabel('Time (sec)')
        pl.ylabel('Instantaneous HR (BPM)')
        #pl.title('Dynamic Pacing Protocol to mimic HRV')
        pl.legend([HF_model_label + protocol_label])

    # If multiple stimuli in one AP, messes up using offsets, instead record starts of APs and pacing at that start
    offset_array = np.asarray(offset_array)
    offset_array1 = np.insert(offset_array,0,0)
    pcl_start = np.zeros(len(duration))
    for i in range(1,len(duration)):
        z = np.nonzero(offset_array1 < start[i-1])[0][-1]
        pcl_start[i] = pace_array[z]



    ## APD vs time plot ##
    ## ---------------- ##
    pcl_start_sec = [element/1000.0 for element in pcl_start]

    if APD_time_plot == True:
        pl.figure()
        pl.plot(start, duration,'.-')
        pl.xlabel('Time (sec)')
        pl.ylabel('APD (ms)')
        #pl.title('APD varying over time with the HRV protocol')
        #pl.xlim(500,650)
    # Show plots
    pl.show()

    return(pcl_start, duration)
# Main function for testing
def main():
    '''
    pcl_APs, duration =  HRV_return(model = 'ohara-2011', number_runs = 30, HF_protocol = None, HF_model = None, cell_type = 0,  AP_plot = False, protocol_plot = False)
    hf_pcl_APs, hf_duration =  HRV_return(model = 'ohara-2011', number_runs = 30, HF_protocol = 'HF', HF_model = 'Elshrif', cell_type = 0,  AP_plot = False, protocol_plot = False)
    hf_walk_pcl_APs, hf_walk_duration =  HRV_return(model = 'ohara-2011', number_runs = 30, HF_protocol = 'HF_walk', HF_model = 'Elshrif', cell_type = 0,  AP_plot = False, protocol_plot = False)

    a_pcl_APs, a_duration =  HRV_return(model = 'ohara-2011', number_runs = 30, HF_protocol = 'HF', HF_model = None, cell_type = 0,  AP_plot = False, protocol_plot = False)
    b_pcl_APs, b_duration =  HRV_return(model = 'ohara-2011', number_runs = 30, HF_protocol = None, HF_model = 'Elshrif', cell_type = 0,  AP_plot = False, protocol_plot = False)

    pl.plot(pcl_APs[len(pcl_APs)/2 :], duration[len(pcl_APs)/2:], '.')
    pl.plot(hf_pcl_APs[len(hf_pcl_APs)/2 :], hf_duration[len(hf_pcl_APs)/2:], '.')
    pl.plot(hf_walk_pcl_APs[len(hf_walk_pcl_APs)/2 :], hf_walk_duration[len(hf_walk_pcl_APs)/2:], '.')

    pl.plot(a_pcl_APs[len(a_pcl_APs)/2 :], a_duration[len(a_pcl_APs)/2:], '.')
    pl.plot(b_pcl_APs[len(b_pcl_APs)/2 :], b_duration[len(b_pcl_APs)/2:], '.')

    pl.legend(['Normal model, healthy protocol', 'HF model, HF protocol at rest', 'HF model, HF protocol 6-min walk test','Normal model, HF protocol', 'HF model, healthy protocol'])
    '''
    pcl_APs, duration =  HRV_return(model = 'ohara-2011', number_runs = 50, HF_protocol = 'HF', HF_model = 'Elshrif', cell_type = 0,  AP_plot = True, protocol_plot = False, min_PCL = 350, max_PCL =480, APD_time_plot = True)
    pl.plot(pcl_APs[len(pcl_APs)/2 :], duration[len(pcl_APs)/2:], '.')

    pl.show()
if __name__ == "__main__":
    main()
