#!/usr/bin/env python
import myokit
import matplotlib.pyplot as pl
from manual_APD import ap_duration
from apd_dynamic import apd_dynamic
import numpy as np
from HF_model import *

## --------------------------------------------- ##
## Functions for different restitution protocols ##
## --------------------------------------------- ##

# Easy to compare with standardised version for all models and cell types

## Dynamic Protocol ##
## ---------------- ##

# Plots PCL vs APD (repolarisation % to be specified)
# Specify model and cell type

#models = ['tentusscher-2006', 'grandi-2010', 'ohara-2011', 'ohara-cipa-v1-2017', 'grandi-2010_modified']
#HF_model = ['Gomez','Elshrif','Moreno', 'Lu']
def dynamic_protocol(model, time_per_stage = 30000, stimuli_per_pace = None , max_pcl = 1000, min_pcl = 50, number_stages = 20, repolarisation = 90, cell_type = 1, voltage_plot = False, HF_model = None, plot_from_function = True):
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

    p = myokit.Protocol()
    m = myokit.load_model('{}.mmt'.format(model))
    models_HF = {'Ord_HF_Gomez' : Ord_HF_Gomez, 'Ord_HF_Elshrif' : Ord_HF_Elshrif, 'GPB_HF_Gomez' : GPB_HF_Gomez, 'GPB_HF_Gomez': GPB_HF_Moreno, 'TT_HF_Lu' : TT_HF_Lu}
    if HF_model != None:
        m_str = '{}_HF_{}'.format(HF_label, HF_model)
        m_func = models_HF[m_str]
        hf_m = m_func(cell_type)

    # Empty arrays to fill.
    period = []
    offset_list = []
    offset = 0

    # More points around high pacing, more likely to see graph bifurcate
    pacing_list = np.logspace(np.log10(min_pcl), np.log10(max_pcl), num = number_stages, base = 10.0)[::-1]
    pacing_list = [int(i) for i in pacing_list]

    # Starting at low frequency pacing (1Hz), 30 seconds at each pacing and then moving towards 230ms
    for pacing in pacing_list:
        beats_per_pace = (time_per_stage)/pacing
        if stimuli_per_pace != None:
            beats_per_pace = stimuli_per_pace
            no_runs = 50*np.sum(pacing_list)
        else:
            no_runs = time_per_stage*len(pacing_list)
        p.schedule(1, start = offset, duration =0.5, period = pacing, multiplier = beats_per_pace)
        offset_list.append(offset)
        # Next set of pacing events to be scheduled by this offset (beats per pace * pace)
        offset += beats_per_pace*pacing
        period.append(pacing)

    # Set up simulation using this scheduled protocol
    s = myokit.Simulation(m, p)
    s.set_constant(label, cell_type)

    ## HF model as well ##
    if HF_model != None:
        hf_s = myokit.Simulation(hf_m, p)
        hf_s.set_constant(label, cell_type)
        hf_d = hf_s.run(offset, log = ['membrane.V','engine.time'])
        hf_start, hf_duration, hf_thresh = apd_dynamic(hf_d, p, no_runs, repolarisation = repolarisation)

        hf_final_apd = np.zeros(len(offset_list))
        hf_final_apd2 = np.zeros(len(offset_list))
        hf_final_apd3 = np.zeros(len(offset_list))
        hf_final_apd4 = np.zeros(len(offset_list))

    # Run the simulation with final offset value, equal to time passed for whole protocol
    d = s.run(offset, log = ['membrane.V','engine.time'])

    # Including log of paces, as pacing is not constant, changes threshold for each pace
    logint = p.create_log_for_interval(0, p.characteristic_time(), for_drawing = True)

    # Use ap_duration function to calculate start times and durations
    #start, duration, thresh = ap_duration(d, 30000*len(pacing_list), repolarisation = 95, log_for_interval = logint)
    start, duration, thresh = apd_dynamic(d, p, no_runs, repolarisation = repolarisation)

    # Numpy array to contain final APD for each pacing cycle (and 3 previous APs in case of alternans)
    final_apd = np.zeros(len(offset_list))
    final_apd2 = np.zeros(len(offset_list))
    final_apd3 = np.zeros(len(offset_list))
    final_apd4 = np.zeros(len(offset_list))

    # First offset equal to zero, so remove first entry from the list
    offset_list = offset_list[1:]

    # If user wants a plot to check individual APs
    if voltage_plot == True:
        # Plot time vs membrane potential graph. Check right APs are being indexed
        pl.figure()
        pl.plot(d['engine.time'],d['membrane.V'])
        pl.xlabel('Time (ms)')
        pl.ylabel('Membrane Potential (mV)')
        for value in offset_list:
            pl.axvline(x= value, color = 'red', ls = 'dotted', ymin = 0.05, ymax = 0.96)


    # Fill arrays with APDs for last few APs at each PCL
    for i in range(len(offset_list)):
        # Final peak of pacing cycle = peak before the first of a new pacing cycle

        # Index_start = The index of the start of pacing cycle in start array
        index_start = np.nonzero(start >= offset_list[i])[0][0]

        pl.arrow(start[index_start-2], thresh[index_start-2], duration[index_start-2], 0, head_width=2, head_length=duration[index_start-2]/3, length_includes_head=True)
        #pl.text(start + duration_ap/3, -90, str(int(duration_ap)) + ' ms')

        # index_1 -1 to get peak at the end of the previous cycle
        final_apd[i] = duration[index_start -2]
        final_apd2[i] = duration[index_start - 3]
        final_apd3[i] = duration[index_start - 4]
        final_apd4[i] = duration[index_start - 5]
        if HF_model != None:
            hf_index_start = np.nonzero(hf_start >= offset_list[i])[0][0]
            hf_final_apd[i] = hf_duration[hf_index_start -2]
            hf_final_apd2[i] = hf_duration[hf_index_start - 3]
            hf_final_apd3[i] = hf_duration[hf_index_start - 4]
            hf_final_apd4[i] = hf_duration[hf_index_start - 5]

    # The final peak doesn't have an offset after, so can be indexed by the final duration recorded (take a few for alternans)
    final_apd[-1] = duration[-2]
    final_apd2[-1] = duration[-3]
    final_apd3[-1] = duration[-4]
    final_apd4[-1] = duration[-5]

    if HF_model != None:
        hf_final_apd[-1] = hf_duration[-2]
        hf_final_apd2[-1] = hf_duration[-3]
        hf_final_apd3[-1] = hf_duration[-4]
        hf_final_apd4[-1] = hf_duration[-5]

    # Plot the restitution curve
    pl.figure()
    pl.plot(pacing_list, final_apd, '.', c = 'b')
    pl.plot(pacing_list,final_apd2,'.', c = 'b')
    pl.plot(pacing_list,final_apd3,'.', c = 'b')
    pl.plot(pacing_list,final_apd4,'.', c = 'b')

    if HF_model != None:
        hf_final_apd[-1] = hf_duration[-2]
        hf_final_apd2[-1] = hf_duration[-3]
        hf_final_apd3[-1] = hf_duration[-4]
        hf_final_apd4[-1] = hf_duration[-5]
        pl.plot(pacing_list, hf_final_apd, '.', c = 'C1')
        pl.plot(pacing_list, hf_final_apd2,'.', c = 'C1')
        pl.plot(pacing_list, hf_final_apd3,'.', c = 'C1')
        pl.plot(pacing_list, hf_final_apd4,'.', c = 'C1')
        pl.legend(['Healthy', 'HF'])

    pl.xlabel('PCL (ms)')
    pl.ylabel('APD {} (ms)'.format(repolarisation))

    pl.title('{} {} Cells Dynamic Protocol Restitution Curve'.format(name,cell_types[cell_type]))
    pl.show()

## S1- S2 Protocol ##
## --------------- ##

#models = ['tentusscher-2006', 'grandi-2010', 'ohara-2011', 'ohara-cipa-v1-2017.mmt']
def s1s2_protocol(model, number_S1, PCL_S1 = 1000, pre_pacing = 200, min_di = 10, max_di = 1000, number_di = 30, repolarisation = 90, cell_type = 1, log_scale = False, HF_model = None, plot_from_function = True):
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

    p = myokit.Protocol()
    m = myokit.load_model('{}.mmt'.format(model))
    models_HF = {'Ord_HF_Gomez' : Ord_HF_Gomez, 'Ord_HF_Elshrif' : Ord_HF_Elshrif, 'GPB_HF_Gomez' : GPB_HF_Gomez, 'GPB_HF_Gomez': GPB_HF_Moreno, 'TT_HF_Lu' : TT_HF_Lu, 'Ordcipa_HF_Gomez' : Ordcipa_HF_Gomez, 'Ordcipa_HF_Elshrif' : Ordcipa_HF_Elshrif}
    if HF_model != None:
        m_str = '{}_HF_{}'.format(HF_label, HF_model)
        m_func = models_HF[m_str]
        m = m_func(cell_type)

    # Set indefinitely recurring event of constant bcl
    p = myokit.pacing.blocktrain(PCL_S1, 0.5, offset=0, level=1.0, limit=0)
    s = myokit.Simulation(m, p)
    s.set_constant(label, cell_type)

    # Pre-pace with these conditions
    s.pre(pre_pacing*PCL_S1)
    # Set number of S1 beats to record after pre-pacing
    p = myokit.pacing.blocktrain(PCL_S1, 0.5, offset=0, level=1.0, limit = number_S1)
    s.set_protocol(p)
    d = s.run(number_S1*PCL_S1)

    # Calculate end of AP for final S1 beat (APD 90). Can add DI to this value for S2 start time
    start, duration, thresh = ap_duration(d, number_S1, repolarisation = 95)
    end_final_s1 = start[-1] + duration[-1]

    # Reset
    s.reset()

    di_list = np.zeros(number_di)
    apd_list = np.zeros(number_di)
    i = 0
    # Equally space points on base 10 log scale
    if log_scale == True:
        interval = np.logspace(int(np.log10(min_di)), int(np.log10(max_di)), num = number_di, base = 10.0)
    else:
        interval = np.linspace(min_di, max_di, num = number_di)
    for di in interval:

        # Finite number of S1 beats
        p = myokit.pacing.blocktrain(PCL_S1, 0.5, offset=0, level=1.0, limit=number_S1)
        # Schedule S2 beat for the end of last AP (using APD 90) + given DI
        p.schedule(1, end_final_s1 + di, 0.5, PCL_S1, 1)

        # Set protocol, run and calculate APDs
        s.set_protocol(p)
        d = s.run((number_S1 + 1)*PCL_S1 + di)
        start, duration, thresh = ap_duration(d, number_S1 + 2, repolarisation = repolarisation)

        # Storing apd and DI for this pacing length
        apd_list[i] = (duration[-1])
        di_list[i] = (di)
        s.reset()
        i += 1
    grad_list = np.zeros(len(di_list) -1)
    for i in range(1,len(di_list)):
        grad_list[i-1] =  (apd_list[i]-apd_list[i-1])/(di_list[i] - di_list[i-1])
    max_grad = max(grad_list)

    if plot_from_function == True:
        pl.figure()
        # Plot line once iterated over all DI values
        pl.plot(di_list,apd_list)
        pl.plot(di_list, apd_list, 'x', label = '_nolegend_')

        # Plot S1S2 protocol restitution curve
        pl.xlabel('Diastole interval (ms)')
        pl.ylabel('APD {}(ms)'.format(repolarisation))
        pl.title('{} {} Cell, {} ms PCL S1-S2 Protocol Restitution Curve'.format(name, cell_types[cell_type], PCL_S1))
        if log_scale == True:
            pl.xscale('log')
            pl.xlim(min_di, max_di)
        pl.show()

    return(di_list, apd_list, max_grad)

## Multistability Test Protocol. S1-CI-S2. ##
## --------------------------------------- ##

def s1cis2(model, time_per_stage, max_pcl = 1000, min_pcl = 100, number_stages = 20, repolarisation = 90, cell_type = 1, voltage_plot = False, CIs = [120,150]):
    # S1-CI-C2. User should have defined number of bifurcations plots to produce, by number of CIs given.
    number_CI = len(CIs)
    CIs = np.asarray(CIs)
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

    # Log spaced : more points around high pacing, more likely to see graph bifurcate
    #pacing_list = np.logspace(np.log10(min_pcl), np.log10(max_pcl), num = number_stages, base = 10.0)[::-1]
    # Linearly spaced points
    pacing_list = np.linspace(min_pcl, max_pcl, num = number_stages)[::-1]
    # Only integer values used
    pacing_list = [int(i) for i in pacing_list]

    for ci in CIs:
        # Empty arrays to fill.
        period = []
        offset_list = []
        offset = 0
        # Blank protocol to add to
        p = myokit.Protocol()

        # Starting at low frequency pacing (1Hz), 30 seconds at each pacing and then moving towards 230ms
        for pacing in pacing_list:
            beats_per_pace = (time_per_stage)/pacing
            p.schedule(1, start = offset, duration =0.5, period = pacing, multiplier = beats_per_pace)
            offset_list.append(offset)
            # Next set of pacing events to be scheduled, CI from last stimulus of S1 beat to first stimulus of S2 beat
            offset += (beats_per_pace-1)*pacing + ci
            period.append(pacing)

        # Set up simulation using this scheduled protocol
        s = myokit.Simulation(m, p)
        s.set_constant(label, cell_type)

        # Run the simulation with final offset value, equal to time passed for whole protocol
        d = s.run(offset, log = ['membrane.V','engine.time'])

        # Including log of paces, as pacing is not constant, changes threshold for each pace
        logint = p.create_log_for_interval(0, p.characteristic_time(), for_drawing = True)

        # Use ap_duration function to calculate start times and durations
        #start, duration, thresh = ap_duration(d, 30000*len(pacing_list), repolarisation = 95, log_for_interval = logint)
        start, duration, thresh = apd_dynamic(d, p, time_per_stage*len(pacing_list), repolarisation = repolarisation)

        # Numpy array to contain final APD for each pacing cycle (and 3 previous APs in case of alternans)
        final_apd = np.zeros(len(offset_list))
        final_apd2 = np.zeros(len(offset_list))
        final_apd3 = np.zeros(len(offset_list))
        final_apd4 = np.zeros(len(offset_list))

        # First offset equal to zero, so remove first entry from the list
        offset_list = offset_list[1:]

        # If user wants a plot to check individual APs
        if voltage_plot == True:
            # Plot time vs membrane potential graph. Check right APs are being indexed
            pl.figure()
            pl.plot(d['engine.time'],d['membrane.V'])
            pl.xlabel('Time (ms)')
            pl.ylabel('Membrane Potential (mV)')
            for value in offset_list:
                pl.axvline(x= value, color = 'red', ls = 'dotted', ymin = 0.05, ymax = 0.96)


        # Fill arrays with APDs for last few APs at each PCL
        for i in range(len(offset_list)):
            # Final peak of pacing cycle = peak before the first of a new pacing cycle

            # Index_start = The index of the start of pacing cycle in start array
            print start
            print offset_list
            index_start = np.nonzero(start >= offset_list[i])[0][0]

            pl.arrow(start[index_start-2], thresh[index_start-2], duration[index_start-2], 0, head_width=2, head_length=duration[index_start-2]/3, length_includes_head=True)
            #pl.text(start + duration_ap/3, -90, str(int(duration_ap)) + ' ms')

            # index_1 -1 to get peak at the end of the previous cycle
            final_apd[i] = duration[index_start -2]
            final_apd2[i] = duration[index_start - 3]
            final_apd3[i] = duration[index_start - 4]
            final_apd4[i] = duration[index_start - 5]

        # The final peak doesn't have a peak after it, so can be indexed by the final duration recorded
        final_apd[-1] = duration[-2]
        # If the graph has a long-short pattern, take previous APDs and DIs as well
        final_apd2[-1] = duration[-3]
        final_apd3[-1] = duration[-4]
        final_apd4[-1] = duration[-5]

        # Plot the restitution curve
        pl.figure()
        pl.plot(pacing_list, final_apd, 'x', c = 'b')
        pl.plot(pacing_list,final_apd2,'x', c = 'b')
        pl.plot(pacing_list,final_apd3,'x', c = 'b')
        pl.plot(pacing_list,final_apd4,'x', c = 'b')
        pl.xlabel('PCL (ms)')
        pl.ylabel('APD {} (ms)'.format(repolarisation))
        pl.title('{} {} Cells Multistability Test Protocol Bifurcation Diagram. CI = {} ms'.format(name,cell_types[cell_type], ci))

    pl.show()

#s1s2_protocol('ohara-2011', number_S1 = 10, PCL_S1 = 1000, pre_pacing = 200, min_di = 10, max_di = 1000, number_di = 50, repolarisation = 90, cell_type = 1, log_scale = False, HF_model = 'Elshrif')
