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

offset_list = []
pcl = 600
p = myokit.pacing.blocktrain(pcl, 0.5, offset=20, level=1.0, limit=0)
s = myokit.Simulation(m,p)
s.pre(200*pcl)
p = myokit.Protocol()
offset = 20

interval = [300,310,350,390,430,470,510,550,590,670,750,830,950,1000,1100,1250,1400,1500, 1600][::-1]
for lag in interval:
    beats_per_pace = 20
    p.schedule(1,start = offset, duration =0.5, period = pcl, multiplier = beats_per_pace)
    offset_list.append(offset)
    # Next set of lag events to be scheduled by this offset (beats per pace * pace)
    offset += (beats_per_pace-1)*pcl + lag



# Set up simulation using this scheduled protocol
s.set_protocol(p)
s.set_constant('type.epi', cell_types[cell_type])
# Run the simulation with final offset value, equal to time passed for whole protocol
d = s.run(offset, log = ['membrane.V','engine.time'])
print 'finished running'

# Use ap_duration function to calculate start times and durations
start, duration, thresh = ap_duration(d, beats_per_pace*len(interval))
print 'apd calculated'

# First offset equal to zero, so remove first entry from the list
offset_list = offset_list[1:]

# Numpy array to contain final APD for each lag cycle
final_apd = np.zeros(len(offset_list) + 1)

index_start_list = []
# Fill final_apd array
for i in range(len(offset_list)):
    # Final peak of lag cycle = peak before the first of a new lag cycle

    # Index_start = The index of the start of lag cycle in start array
    index_start = np.nonzero(start >= offset_list[i])[0][0]
    index_start_list.append(index_start)
    # index_1 -1 to get peak at the end of the previous cycle
    final_apd[i] = duration[index_start]

# The final peak doesn't have a peak after it, so can be indexed by the final duration recorded
final_apd[-1] = duration[-1]


pl.figure()
pl.plot(d['engine.time'],d['membrane.V'])
for i, in_start in enumerate(index_start_list):
    pl.arrow(start[in_start], -80, final_apd[i], 0, head_width=5, head_length=100,
    length_includes_head=True, color = 'cyan')

# Plot the restitution curve
pl.figure()
pl.plot(interval, final_apd, 'x', c = 'b')
pl.xlabel('Extrastimulus interval (ms)')
pl.ylabel('APD 90 (ms)')
pl.ylim(200,350)
pl.title('Dynamic Restitution Curve- Grandi (2010)')
pl.show()
