#!/usr/bin/env python
import myokit
import matplotlib.pyplot as pl
from manual_APD import ap_duration
import numpy as np
from HF_model import Ord_HF_Gomez

## Protocol same for all models and cell types, to compare ##
## Use O'hara dynamic protocol (with PCL not DI) and ten-Tusscher S1-S2 with more pacing ##

# Using APD 90
percent = 90
bcl = 600
cell_types = {0:'Endocardial', 1: 'Epicardial'}
models = ['tentusscher-2006', 'grandi-2010', 'ohara-2011']
i = 1
pl.figure()
for model in models:
    m = myokit.load_model('{}.mmt'.format(model))
    p = myokit.Protocol()
    s = myokit.Simulation(m,p)
    if model == 'tentusscher-2006':
        label = 'cell.type'
    elif model == 'grandi-2010':
        label = 'type.epi'
    elif model == 'ohara-2011':
        label = 'cell.mode'
    for cell_type in range(0,2):
        s.set_constant(label, cell_type)
        # S1 S2 Protocol. Basic drive length 600ms
        p = myokit.pacing.blocktrain(bcl, 0.5, offset=0, level=1.0, limit=0)
        s.set_protocol(p)
        s.pre(bcl * 200)
        p = myokit.Protocol()
        paces  = 3
        p = myokit.pacing.blocktrain(bcl, 0.5, offset=0, level=1.0, limit=paces)
        d = s.run(paces*bcl)
        start, duration, thresh = ap_duration(d, paces, repolarisation = percent)
        end_final_s1 = start[-1] + duration[-1]
        s.reset()
        pl.subplot(3,2,i)
        di_list = []
        apd_duration = []
        for di in range(10,1000,50):
            p = myokit.pacing.blocktrain(bcl, 0.5, offset=0, level=1.0, limit=paces)
            p.schedule(1, end_final_s1 + di, 0.5, bcl, 1)
            s.set_protocol(p)
            d = s.run((paces + 2)*bcl)
            start, duration, thresh = ap_duration(d, paces, repolarisation = percent)

            # Storing apd and DI for this pacing length
            apd_duration.append(duration[-1])
            di_list.append(di)

            s.reset()

        pl.plot(di_list,apd_duration)
        pl.plot(di_list,apd_duration, 'x', label = '_nolegend_')
        pl.xlabel('Diastole interval (ms)')
        pl.ylabel('APD (ms)')
        pl.text(300, (apd_duration[0]+apd_duration[-1])/2,'{} model, {} cell'.format(model, cell_types[cell_type]))
        i += 1
pl.suptitle('S1-S2 Restitution Protocol, pre-paced for 20 beats at 600ms, 3 S1 beats followed by one S2 beat')
pl.show()
