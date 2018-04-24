#!/usr/bin/env python
import myokit
import matplotlib.pyplot as pl
from manual_APD import ap_duration

### Ten-Tusscher S1-S2 protocol- Figure 5B ###
### -------------------------------------- ###

## 10 S1 beats at 600ms PCL, followed by single S2 beat with varying pacing
# Get the model and protocol, create a simulation
m = myokit.load_model('tentusscher-2006.mmt')

# Pacing for S1 beats is 600ms
pcl = 600
p = myokit.pacing.blocktrain(pcl, 0.5, offset=20, level=1.0, limit=0)
s = myokit.Simulation(m, p, apd_var='membrane.V')

# Set cell type
cell_types = {'Endocardial': 0, 'Epicardial' : 1, 'Mid-myocardial' : 2}
cell_type = 'Epicardial'
s.set_constant('cell.type', cell_types[cell_type])

#Pre pace 9 beats (10 s1 beats in total)
s.pre(pcl * 9)

# Setting a step size to increase pacing for S2 interval by
step_size = 40

# Final S1 beat and S2 beat, with APD 50 and APD 90
# PCL = APD +DI
pl.figure()
# For loop for plotting APD 50% and APD 90%
for percent in range(50,91,40):
    di = []
    apd_duration= []
    pacing = 342
    while pacing <1000:
        #print pacing

        # Final S1 beat
        p = myokit.pacing.blocktrain(pcl, 0.5, offset=0, level=1.0, limit=1)

        # Introducing S2 upstroke. Offset by pacing interval
        p.schedule(1, pacing,0.5, pcl, 1)
        s = myokit.Simulation(m, p,apd_var='membrane.V')

        # Run actual simulation to calculate APD for
        paces = 5
        d = s.run(paces*pcl)

        # Run using function.
        start, duration = ap_duration(d, paces, repolarisation = percent)

        # Storing apd and DI for this pacing length
        apd_duration.append(duration[1])
        di.append(pacing - duration[0])

         #Step size of pacing interval increase
        pacing += step_size

        # Reset simulation settings to pre-pacing
        s.reset()

    pl.plot(di, apd_duration)
    pl.plot(di, apd_duration, 'x',label='_nolegend_')

# Plot S1S2 protocol restitution curve
pl.xlabel('Diastole interval (ms)')
pl.ylabel('APD (ms)')
pl.title('Ten-Tusscher (2006) {} Cell 10 x S1, 1 x S2 Protocol Restitution Curve (Fig 5B)'.format(cell_type))
pl.legend(['APD 50','APD 90'])
pl.xlim(0,600)
pl.show()

# Parameters matching 2nd row table 2 ten-Tusscher model 2006 (should be slope of 1.1)
