#!/usr/bin/env python

import myokit
import matplotlib.pyplot as pl
from manual_APD import ap_duration

### Ten-Tusscher Dynamic Restitution Protocol- Figure 5F/J ###
### ------------------------------------------------------ ###
# Starting from low pacing and increasing frequency step wise
# Measuring the final AP for each pacing

# Get the model
m = myokit.load_model('tentusscher-2006.mmt')

step_size = -30
apd_duration = []
alternate_apd = []
#di = []
period = []
pcl_list  = []

pcl = 800
p = myokit.pacing.blocktrain(pcl, 0.5, offset=0, level=1.0, limit=400)
s = myokit.Simulation(m, p,apd_var='membrane.V')
s.pre(pcl * 400)
s.reset()
offset = 0
while pcl > 100:
    print pcl
    #Initialising protocol with new pacing
    '''
    p.schedule(1, offset, 0.5, pcl, 50 )

    offset += 50*pcl
    pcl += step_size
    '''
    p = myokit.pacing.blocktrain(pcl, 0.5, offset=0, level=1.0, limit=0)
    s = myokit.Simulation(m, p,apd_var='membrane.V')

    #Pre pacing for 48 beats, running 10 beats
    s.pre(pcl * 40)
    thresh = 0.9 * s.state()[m.get('membrane.V').indice()]
    d = s.run(10*pcl)
    start, duration = ap_duration(d, 10, thresh)

    apd_duration.append(duration[-1])
    alternate_apd.append(duration[-2])
    period.append(pcl)
    pcl += step_size
    s.reset()

# Plotting APD vs Period for dynamic protocol
pl.figure()
pl.xlabel('Period (ms)')
pl.ylabel('APD 90 (ms)')
pl.title('Ten-Tusscher (2006)- Dynamic restitution curve (Fig. 5J)')
pl.xlim(0,800)

pl.plot(period, apd_duration, 'x')
pl.plot(period, alternate_apd, 'x')
pl.show()
