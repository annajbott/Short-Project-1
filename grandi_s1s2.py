#!/usr/bin/env python
import myokit
import matplotlib.pyplot as pl
from manual_APD import ap_duration
model = 'grandi-2010'
s1s2_protocol(model, number_S1, PCL_S1 = 1000, pre_pacing = 200, min_di = 10, max_di = 1000, number_di = 30, repolarisation = 90, cell_type = 0, log_scale = False, HF_model = None, plot_from_function = True)


# run 
