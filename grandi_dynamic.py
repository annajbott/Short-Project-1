#!/usr/bin/env python
import myokit
import matplotlib.pyplot as pl
from manual_APD import ap_duration
import numpy as np

# Grandi Dynamic Restitution ##

'''
APD was measured as the interval between AP upstroke and 90% repolarization level
(APD 90). APDrestitution was assessed with a dynamic protocol: after establishing a
steady-state at a basic cycle length of 600 ms, extrastimuli were applied at 20
cycle intervals, reducing the coupling interval in decrements of 40, 20 and
10 ms from the drive cycle length. Stimuli were also delivered at coupling
intervals larger than the steady-state drive cycle length as in Morgan et al.
'''

m = myokit.load_model('grandi-2010.mmt')
cell_types = {'Endocardial': 0, 'Epicardial' : 1}
cell_type = 'Endocardial'

# Empty arrays to fill.
period = []
di = []
offset_list = []

p = myokit.Protocol()
offset = 0
# More points around high pacing, more likely to see graph bifurcate
pacing_list = [230,250,270,290,310,335,350,390,430,470,510,550,590,670,750,830,950,1000, 1400, 1600][::-1]
# Starting at low frequency pacing (1Hz), 30 seconds at each pacing and then moving towards 230ms
for pacing in pacing_list:
    beats_per_pace = 100
    p.schedule(1,start = offset, duration =0.5, period = pacing, multiplier = beats_per_pace)
    offset_list.append(offset)
    # Next set of pacing events to be scheduled by this offset (beats per pace * pace)
    offset += beats_per_pace*pacing
    period.append(pacing)


# Set up simulation using this scheduled protocol
s = myokit.Simulation(m, p)
s.set_constant('type.epi', cell_types[cell_type])
# Run the simulation with final offset value, equal to time passed for whole protocol
d = s.run(offset, log = ['membrane.V','engine.time'])
print 'finished running'

# Use ap_duration function to calculate start times and durations
start, duration, thresh = ap_duration(d, beats_per_pace*len(pacing_list))
print 'apd calculated'

# First offset equal to zero, so remove first entry from the list
offset_list = offset_list[1:]
print offset_list
print start

# Numpy array to contain final APD for each pacing cycle
final_apd = np.zeros(len(offset_list) + 1)
final_apd2 = np.zeros(len(offset_list) + 1)
final_apd3 = np.zeros(len(offset_list) + 1)
final_apd4 = np.zeros(len(offset_list) + 1)
di2 = []
di3 = []
di4 = []

# Fill final_apd array
for i in range(len(offset_list)):
    # Final peak of pacing cycle = peak before the first of a new pacing cycle

    # Index_start = The index of the start of pacing cycle in start array
    index_start = np.nonzero(start >= offset_list[i])[0][0]

    # index_1 -1 to get peak at the end of the previous cycle
    final_apd[i] = duration[index_start-1]
    final_apd2[i] = duration[index_start - 2]
    final_apd3[i] = duration[index_start - 3]
    final_apd4[i] = duration[index_start - 4]
    di.append(period[i]-duration[index_start -2])
    di2.append(period[i]-duration[index_start -3])
    di3.append(period[i]-duration[index_start -4])
    di4.append(period[i]-duration[index_start -5])

# The final peak doesn't have a peak after it, so can be indexed by the final duration recorded
final_apd[-1] = duration[-1]
di.append(period[-1] - duration[-2])

# If the graph has a long-short pattern, take previous APD and DI as well
final_apd2[-1] = duration[-2]
di2.append(period[-1] - duration[-3])
final_apd3[-1] = duration[-2]
di3.append(period[-1] - duration[-3])
final_apd4[-1] = duration[-2]
di4.append(period[-1] - duration[-3])


pl.figure()
pl.plot(d['engine.time'],d['membrane.V'])

# Plot the restitution curve
pl.figure()
pl.plot(pacing_list, final_apd, 'x', c = 'b')
#pl.plot(pacing_list,final_apd2,'x', c = 'b')
#pl.plot(pacing_list,final_apd3,'x', c = 'b')
#pl.plot(pacing_list,final_apd4,'x', c = 'b')
#pl.xlabel('PCL (ms)')
#pl.plot(di, final_apd, 'x', c = 'b')
#pl.plot(di2,final_apd2,'x', c = 'b')
#pl.plot(di3,final_apd3,'x', c = 'b')
#pl.plot(di4,final_apd4,'x', c = 'b')
#pl.xlabel('DI (ms)')

pl.ylabel('APD 90 (ms)')
pl.title('Dynamic Restitution Curve- Grandi (2010)')
pl.show()
