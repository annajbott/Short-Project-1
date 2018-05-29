#!/usr/bin/env python
import myokit
import matplotlib.pyplot as pl
from manual_APD import ap_duration
import numpy as np
from HF_model import Ord_HF_Gomez
from restitution_protocols import *

## Protocol same for all models and cell types, to compare ##
## Use O'hara dynamic protocol (with PCL not DI) and ten-Tusscher S1-S2 with more pacing ##

cell_types = {0:'Endocardial', 1: 'Epicardial'}
models = ['tentusscher-2006', 'grandi-2010_modified', 'ohara-2011', 'ohara-cipa-v1-2017']
i = 1
pl.figure()
for model in models:
    m = myokit.load_model('{}.mmt'.format(model))
    p = myokit.Protocol()
    s = myokit.Simulation(m,p)
    if model == 'tentusscher-2006':
        label = 'cell.type'
        HF_model = 'Lu'
    elif model == 'grandi-2010':
        label = 'type.epi'
        HF_model = 'Gomez'
    elif model == 'grandi-2010_modified':
        label = 'type.epi'
        HF_model = 'Gomez'
    elif model == 'ohara-2011':
        label = 'cell.mode'
        HF_model = 'Gomez'
    elif model == 'ohara-cipa-v1-2017':
        label = 'cell.celltype'
        HF_model = 'Gomez'
    # Not mid-myocardial cells
    for cell_type in range(0,2):
        # Each model on the same line
        pl.subplot(4,2,i)

        # Run S1S2 protocol using function, turn internal function plot off so can plot using subplots here instead
        di_list, apd_list, max_grad = s1s2_protocol(model, number_S1 = 5, PCL_S1 = 1000, pre_pacing = 50, min_di = 10, max_di = 1000, number_di = 50, repolarisation = 90, cell_type = cell_type, log_scale = False, plot_from_function = False)

        #pl.plot(di_list,apd_list)
        pl.plot(di_list,apd_list, 'k.-')
        pl.text(300, (apd_list[0]+apd_list[-1])/2,'{} model, {} cell'.format(model, cell_types[cell_type]))

        # Plot HF S1S2 curve
        di_list, apd_list, max_grad = s1s2_protocol(model, number_S1 = 5, PCL_S1 = 1000, pre_pacing = 50, min_di = 10, max_di = 1000, number_di = 30, repolarisation = 90, cell_type = cell_type, log_scale = False, HF_model = HF_model, plot_from_function = False)
        pl.plot(di_list,apd_list, 'r.-')

        if model == 'ohara-2011' or model == 'ohara-cipa-v1-2017':
            di_list, apd_list, max_grad = s1s2_protocol(model, number_S1 = 5, PCL_S1 = 1000, pre_pacing = 50, min_di = 10, max_di = 1000, number_di = 30, repolarisation = 90, cell_type = cell_type, log_scale = False, HF_model = 'Elshrif', plot_from_function = False)
            pl.plot(di_list,apd_list, 'c.-')
            pl.legend(['Normal','HF {} '.format(HF_model), 'HF Elshrif'])
        else:
            pl.legend(['Normal','HF {} '.format(HF_model)])

        pl.xlabel('Diastole interval (ms)')
        pl.ylabel('APD (ms)')
        #pl.legend(['Normal','HF {} '.format(HF_model)])
        i += 1

pl.suptitle('S1-S2 Restitution Protocol, pre-paced for 50 beats at 1000ms, 5 S1 beats followed by one S2 beat')
pl.show()
