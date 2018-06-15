#!/usr/bin/env python
import myokit
import matplotlib.pyplot as pl
from manual_APD import ap_duration
import numpy as np
from HF_model import *


## Heart Rate Variability (HRV)- Healthy, HF (rest), HF (submaximal exercise) ##
## -------------------------------------------------------------------------- ##

# Plots for protocol, APs and restitution curve
def HRV_return(model, HF_model = None, HF_protocol = None, number_points_up = 30,number_runs = 50, cell_type = 0, restitution_curve = False , AP_plot = False, protocol_plot = True, max_PCL = None, min_PCL = None, average_init = None, APD_time_plot = False):

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
    models_HF_col = {'Ord_HF_Gomez' : 'orange', 'Ord_HF_Elshrif' : 'limegreen', 'GPB_HF_Gomez' : 'orange', 'GPB_HF_Moreno': 'm', 'TT_HF_Lu' : 'gold', 'Ordcipa_HF_Elshrif' : 'limegreen', 'Ordcipa_HF_Gomez': 'orange'}
    if HF_model != None:
        m_str = '{}_HF_{}'.format(HF_label, HF_model)
        m_func = models_HF[m_str]
        m = m_func(cell_type)
        HF_model_label = 'HF ' + HF_model


    # Protocol
    p = myokit.Protocol()

    # Using sin wave to simulate step protocol
    points = np.linspace(0,2*np.pi, num  = number_points_up, endpoint = False)
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

    # If min and max pcl specified, override protocol options, but only if average not specified
    if max_PCL != None and min_PCL != None and average_init == None:
        max_hr = 60000.0/max_PCL
        min_hr = 60000.0/min_PCL
        average = np.round((max_hr + min_hr)/2.0,0)
        amplitude = np.round((max_hr - min_hr)/2.0,0)

    # If average is specified, take that value, but use amplitudes from earlier
    if average_init != None:
        average = 60000.0/average_init
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
    s.pre(np.sum(pacing_list)*20)
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
    start = [element/1000.0 for element in start]

    if APD_time_plot == True:
        pl.figure()
        if len(start) > len(duration):
            start = start[0:-1]
        pl.plot(start, duration,'.-')
        pl.xlabel('Time (sec)')
        pl.ylabel('APD (ms)')
        #pl.title('APD varying over time with the HRV protocol')
        #pl.xlim(500,650)


    if restitution_curve == True:
        pl.figure()
        pl.plot(pcl_start[4*len(pcl_start)/5 :], duration[4*len(pcl_start)/5:], '.')
        pl.xlabel('PCL (ms)')
        pl.ylabel('APD 90 (ms)')


    # Show plots
    pl.show()

    return(pcl_start, duration)
# Main function for testing
def main():
    #m = 'ohara-2011'
    m = 'ohara-cipa-v1-2017'
    '''
    HF_m = 'Elshrif'
    pcl_APs, duration =  HRV_return(model = m, number_runs = 30, HF_protocol = None, HF_model = None, cell_type = 0,  AP_plot = False, protocol_plot = False)
    hf_pcl_APs, hf_duration =  HRV_return(model = m, number_runs = 30, HF_protocol = 'HF', HF_model = HF_m, cell_type = 0,  AP_plot = False, protocol_plot = False)
    hf_walk_pcl_APs, hf_walk_duration =  HRV_return(model = m, number_runs = 30, HF_protocol = 'HF_walk', HF_model = HF_m, cell_type = 0,  AP_plot = False, protocol_plot = False)

    a_pcl_APs, a_duration =  HRV_return(model = m, number_runs = 30, HF_protocol = 'HF', HF_model = None, cell_type = 0,  AP_plot = False, protocol_plot = False)
    b_pcl_APs, b_duration =  HRV_return(model = m, number_runs = 30, HF_protocol = None, HF_model = HF_m, cell_type = 0,  AP_plot = False, protocol_plot = False)

    pl.plot(pcl_APs[len(pcl_APs)/2 :], duration[len(pcl_APs)/2:], '.', color = 'b')
    pl.plot(hf_pcl_APs[len(hf_pcl_APs)/2 :], hf_duration[len(hf_pcl_APs)/2:], '.', color = 'limegreen')
    pl.plot(hf_walk_pcl_APs[len(hf_walk_pcl_APs)/2 :], hf_walk_duration[len(hf_walk_pcl_APs)/2:], '.', color = 'darkgreen')

    pl.plot(a_pcl_APs[len(a_pcl_APs)/2 :], a_duration[len(a_pcl_APs)/2:], '.', color = 'midnightblue')
    pl.plot(b_pcl_APs[len(b_pcl_APs)/2 :], b_duration[len(b_pcl_APs)/2:], '.', color = 'turquoise')
    pl.legend(['Normal model, healthy protocol', 'HF model, HF rest protocol', 'HF model, HF exercise protocol','Normal model, HF rest protocol', 'HF model, healthy protocol'])
    pl.ylabel('APD 90 (ms)')
    pl.xlabel('PCL (ms)')


    HF_m = 'Gomez'
    pcl_APs, duration =  HRV_return(model = m, number_runs = 30, HF_protocol = None, HF_model = None, cell_type = 0,  AP_plot = False, protocol_plot = False)
    hf_pcl_APs, hf_duration =  HRV_return(model = m, number_runs = 30, HF_protocol = 'HF', HF_model = HF_m, cell_type = 0,  AP_plot = False, protocol_plot = False)
    hf_walk_pcl_APs, hf_walk_duration =  HRV_return(model = m, number_runs = 30, HF_protocol = 'HF_walk', HF_model = HF_m, cell_type = 0,  AP_plot = False, protocol_plot = False)

    a_pcl_APs, a_duration =  HRV_return(model = m, number_runs = 30, HF_protocol = 'HF', HF_model = None, cell_type = 0,  AP_plot = False, protocol_plot = False)
    b_pcl_APs, b_duration =  HRV_return(model = m, number_runs = 30, HF_protocol = None, HF_model = HF_m, cell_type = 0,  AP_plot = False, protocol_plot = False)

    pl.plot(pcl_APs[len(pcl_APs)/2 :], duration[len(pcl_APs)/2:], '.', color = 'b')
    pl.plot(hf_pcl_APs[len(hf_pcl_APs)/2 :], hf_duration[len(hf_pcl_APs)/2:], '.', color = 'orange')
    pl.plot(hf_walk_pcl_APs[len(hf_walk_pcl_APs)/2 :], hf_walk_duration[len(hf_walk_pcl_APs)/2:], '.', color = 'firebrick')

    pl.plot(a_pcl_APs[len(a_pcl_APs)/2 :], a_duration[len(a_pcl_APs)/2:], '.', color = 'midnightblue')
    pl.plot(b_pcl_APs[len(b_pcl_APs)/2 :], b_duration[len(b_pcl_APs)/2:], '.', color = 'darksalmon')
    pl.legend(['Normal model, healthy protocol', 'HF model, HF rest protocol', 'HF model, HF exercise protocol','Normal model, HF rest protocol', 'HF model, healthy protocol'])
    pl.ylabel('APD 90 (ms)')
    pl.xlabel('PCL (ms)')
    '''

    ## Doing cycles in alternan and complex region
    models_ya = ['ohara-2011','ohara-cipa-v1-2017']
    HF_ya = ['Gomez','Elshrif']


    m = models_ya[1]
    HF = HF_ya[1]
    # HF- 1:1 ends 408 (ord-cipa)
    pcl_APs, duration =  HRV_return(model = m, number_runs = 50, HF_protocol = 'HF', HF_model = HF, cell_type = 0, restitution_curve = True,  AP_plot = True, protocol_plot = False, average_init = 408, APD_time_plot = True)

    # No HF- 1:1 end 287 (ord cipa)
    pcl_APs, duration =  HRV_return(model = m, number_runs = 50, HF_protocol = None, HF_model = None, cell_type = 0, restitution_curve = True,  AP_plot = True, protocol_plot = False, average_init = 287, APD_time_plot = True)


    #pl.plot(pcl_APs[2*len(pcl_APs)/3 :], duration[2*len(pcl_APs)/3:], '.')

    #pcl_APs, duration =  HRV_return(model = m, number_points_up = 50, number_runs = 50, HF_protocol = None, HF_model = None, cell_type = 0,  AP_plot = True, protocol_plot = False, min_PCL = 100, max_PCL =600, APD_time_plot = True)
    #pl.plot(pcl_APs[len(pcl_APs)/2 :], duration[len(pcl_APs)/2:], '.')



    pl.show()
if __name__ == "__main__":
    main()
