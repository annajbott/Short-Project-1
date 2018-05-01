#!/usr/bin/env python
import myokit
import matplotlib.pyplot as pl
from manual_APD import ap_duration
import numpy as np

## O'hara Restitution
## S1S2 steady state S1 pacing at 1000ms, single

pcl = 1000
step_size = 50
m = myokit.load_model('ohara-2011.mmt')

p = myokit.pacing.blocktrain(pcl, 0.5, offset=0, level=1.0, limit=0)
s = myokit.Simulation(m, p)

# Set cell type
cell_types = {'Endocardial': 0, 'Epicardial' : 1, 'Mid-myocardial' : 2}
cell_type = 'Endocardial'
s.set_constant('cell.mode', cell_types[cell_type])

#Pre pace to steady state
s.pre(pcl * 300)
apd_90 = []
pacing_list = [280,290,297,298,300,305,310,315,320,330,340,350,380,400,450,550,700,900,1000,3000,5000,7000,9000,10300][::-1]
for percent in [90, 70, 50, 30]:
    di = []
    apd_duration = []
    i = 0
    for pacing in pacing_list:
        p = myokit.Protocol()
        p.schedule(1, 0, 0.5, pcl, 5)
        p.schedule(1, 4*pcl + pacing, 0.5, pcl, 1)
        s.set_protocol(p)
        d = s.run(pacing + 10*pcl)
        start, duration = ap_duration(d, 5, repolarisation = percent)

        apd_duration.append(duration[-1])
        if percent == 90:
            apd_90.append(duration[-2])
        di.append(pacing - apd_90[i])
        i += 1
        s.reset()
    pl.plot(di, apd_duration)
    pl.plot(di, apd_duration, 'x',label='_nolegend_')

# Plot S1S2 protocol restitution curve
pl.xlabel('Diastole interval (ms)')
pl.ylabel('APD (ms)')
pl.title('Ohara {} Cell Steady State S1-S2 Protocol Restitution Curve (Fig 5B)'.format(cell_type))
pl.legend(['APD 90','APD 70','APD 50','APD 30'])
pl.xlim(10,10000)
pl.xscale('log')
pl.show()
'''

### Dynamic Protocol ###
### ---------------- ###

# Empty arrays to fill.
period = []
di = []
offset_list = []

p = myokit.Protocol()
offset = 0
# More points around high pacing, more likely to see graph bifurcate
pacing_list = [230,235,240,245,250,255,260,270,290,310,320,335,350,390,430,470,510,550,590,670,710,750,830,910,950,1000][::-1]
# Starting at low frequency pacing (1Hz), 30 seconds at each pacing and then moving towards 230ms
for pacing in pacing_list:
    beats_per_pace = 30000/pacing
    p.schedule(1,start = offset, duration =0.5, period = pacing, multiplier = beats_per_pace)
    offset_list.append(offset)
    # Next set of pacing events to be scheduled by this offset (beats per pace * pace)
    offset += beats_per_pace*pacing
    period.append(pacing)


# Set up simulation using this scheduled protocol
s = myokit.Simulation(m, p)
s.set_constant('cell.mode', cell_types[cell_type])

# Run the simulation with final offset value, equal to time passed for whole protocol
d = s.run(offset, log = ['membrane.V','engine.time'])

# Use ap_duration function to calculate start times and durations
start, duration = ap_duration(d, 30000*len(pacing_list), repolarisation = 95)

# First offset equal to zero, so remove first entry from the list
offset_list = offset_list[1:]
print offset_list
# Numpy array to contain final APD for each pacing cycle
final_apd = np.zeros(len(offset_list) + 1)
final_apd2 = np.zeros(len(offset_list) + 1)
di2 = []

# Fill final_apd array
for i in range(len(offset_list)):
    # Final peak of pacing cycle = peak before the first of a new pacing cycle

    # Index_start = The index of the start of pacing cycle in start array
    index_start = np.nonzero(np.round(start,0) >= offset_list[i])[0][0]
    print "start", start[index_start], "duration", duration[index_start]
    print "start -1 ", start[index_start -1], "duration", duration[index_start-1]
    print "period", period[i], "di", period[i]-duration[index_start -2]

    # index_1 -1 to get peak at the end of the previous cycle
    final_apd[i] = duration[index_start]
    final_apd2[i] = duration[index_start - 1]
    di.append(period[i]-duration[index_start -1])
    di2.append(period[i]-duration[index_start -2])

# The final peak doesn't have a peak after it, so can be indexed by the final duration recorded
final_apd[-1] = duration[-1]
di.append(period[-1] - duration[-2])

# If the graph has a long-short pattern, take previous APD and DI as well
final_apd2[-1] = duration[-2]
di2.append(period[-1] - duration[-3])

pl.figure()
pl.plot(d['engine.time'],d['membrane.V'])

# Plot the restitution curve
pl.figure()
pl.plot(di, final_apd, 'x', c = 'b')
pl.plot(di2,final_apd2,'x', c = 'b')
pl.xlabel('DI (ms)')
pl.xlim(0,360)
pl.ylim(190,270)
pl.ylabel('APD 95 (ms)')
pl.title('Dynamic Restitution Curve- Endothelial Cells (Ohara 2011)')
pl.show()
'''
