#!/usr/bin/env python
import myokit
import matplotlib.pyplot as pl
import numpy as np
from HF_model import *
from restitution_protocols import *

### Ten-Tusscher Dynamic Restitution Protocol- Figure 5F/J ###
### ------------------------------------------------------ ###
# Starting from low pacing and increasing frequency step wise
# Measuring the final few AP for each pacing

# Get the model
m = myokit.load_model('tentusscher-2006.mmt')
model = 'tentusscher-2006'
#dynamic_protocol(model, number_stages = 80, max_pcl = 1000, min_pcl = 100, voltage_plot = True, cell_type = 0, HF_model = 'Lu', log_scale = False)

#model = 'ohara-2011'

#dynamic_protocol(model, number_stages = 150, time_per_stage = 20000, max_pcl = 1000, min_pcl = 100, HF_model = 'Lu' ,voltage_plot = False, cell_type = 0, log_scale = False)

#model = 'ohara-cipa-v1-2017'
model = 'grandi-2010'
dynamic_protocol(model, number_stages = 150, max_pcl = 1000, min_pcl = 100, HF_model = None ,voltage_plot = False, cell_type = 0, log_scale = False)

#dynamic_protocol(model, number_stages = 60, stimuli_per_pace = 50, min_pcl = 100, voltage_plot = True, cell_type = 0)
'''
di_list, apd_list, max_grad = s1s2_protocol(model, number_S1 = 10, PCL_S1 = 1000, pre_pacing = 200, min_di = 10, max_di = 1000, number_di = 60, repolarisation = 90, cell_type = 0, log_scale = False, HF_model = None, plot_from_function = False)
pl.figure()
pl.plot(di_list, apd_list, '.')

print max_grad

di_list_gomez, apd_list_gomez, hf_g_max_grad = s1s2_protocol(model, number_S1 = 10, PCL_S1 = 1000, pre_pacing = 200, min_di = 10, max_di = 1000, number_di = 60, repolarisation = 90, cell_type = 0, log_scale = False, HF_model = 'Gomez', plot_from_function = False)
di_list_elshrif, apd_list_elshrif, hf_e_max_grad = s1s2_protocol(model, number_S1 = 10, PCL_S1 = 1000, pre_pacing = 200, min_di = 10, max_di = 1000, number_di = 60, repolarisation = 90, cell_type = 0, log_scale = False, HF_model = 'Elshrif', plot_from_function = False)

pl.plot(di_list_gomez, apd_list_gomez, '.')
print hf_g_max_grad
pl.plot(di_list_elshrif, apd_list_elshrif, '.' )
print hf_e_max_grad
pl.legend(['Healthy Ord', 'HF Gomez', 'HF Elshrif'])

pl.ylabel('APD 90 (ms)')
pl.xlabel('DI (ms)')
pl.title('Ord Endocardial Cell S1-S2 Restitution Curve. 1000ms PCL')
pl.show()
'''
