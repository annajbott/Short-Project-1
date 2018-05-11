#!/usr/bin/env python
import myokit
import matplotlib.pyplot as pl
from manual_APD import ap_duration
import numpy as np
from HF_model import Ord_HF_Gomez

## O'hara Restitution
## S1S2 steady state S1 pacing at 1000ms, single
'''
# Using APD 90 for final S1 AP to calculate DI
percent = 90
bcl = 1000
m = myokit.load_model('ohara-2011.mmt')

# HF switch
HF = False
if HF:
    m = Ord_HF_Gomez()

# Set indefinitely recurring event of constant bcl
p = myokit.pacing.blocktrain(bcl, 0.5, offset=0, level=1.0, limit=0)
s = myokit.Simulation(m, p)

# Set cell type
cell_types = {'Endocardial': 0, 'Epicardial' : 1, 'Mid-myocardial' : 2}
cell_type = 'Endocardial'
s.set_constant('cell.mode', cell_types[cell_type])

# Pre-pace with these conditions
s.pre(100*bcl)

# Set number of S1 beats to record after pre-pacing
paces  = 3
p = myokit.pacing.blocktrain(bcl, 0.5, offset=0, level=1.0, limit=paces)
d = s.run(paces*bcl)

# Calculate end of APD for final S1 beat. Can add DI to this value for S2 start time
start, duration, thresh = ap_duration(d, paces, repolarisation = percent)
end_final_s1 = start[-1] + duration[-1]

# Reset so haven't run 3 paces already
s.reset()
pl.figure()

# Have lines on graph for APD 30, 50, 70 and 90
for percent in [90, 70, 50, 30]:
    di_list = []
    apd_duration = []
    # Equally space points on base 10 log scale
    for di in np.logspace(1,4, num = 30, base = 10.0):

        # Finite number of S1 beats
        p = myokit.pacing.blocktrain(bcl, 0.5, offset=0, level=1.0, limit=paces)
        # Schedule S2 beat for the end of last AP (using APD 90) + given DI
        p.schedule(1, end_final_s1 + di, 0.5, bcl, 1)

        # Set protocol, run and calculate APDs
        s.set_protocol(p)
        d = s.run((paces + 1)*bcl + di)
        start, duration, thresh = ap_duration(d, paces + 2, repolarisation = percent)

        # Storing apd and DI for this pacing length
        apd_duration.append(duration[-1])
        di_list.append(di)
        s.reset()

    # Plot line once iterated over all DI values
    pl.plot(di_list,apd_duration)
    pl.plot(di_list,apd_duration, 'x', label = '_nolegend_')

# Plot S1S2 protocol restitution curve
pl.xlabel('Diastole interval (ms)')
pl.ylabel('APD (ms)')
pl.title('Ohara (2011) {} Cell, {} ms PCL S1-S2 Protocol Restitution Curve'.format(cell_type, bcl))
pl.legend(['APD 90','APD 70','APD 50','APD 30'])
pl.xlim(10,10000)
pl.xscale('log')
pl.show()

'''

### Dynamic Protocol ###
### ---------------- ###

m = myokit.load_model('ohara-2011.mmt')
HF = False
if HF:
    m = Ord_HF_Gomez()

# Set cell type
cell_types = {'Endocardial': 0, 'Epicardial' : 1, 'Mid-myocardial' : 2}
cell_type = 'Endocardial'

# Empty arrays to fill.
period = []
di = []
offset_list = []

# Clearing the protocol
p = myokit.Protocol()
offset = 0
# More points around high pacing, more likely to see graph bifurcate
pacing_list = np.logspace(2,3, num = 30, base = 10.0)[::-1]
pacing_list = [int(i) for i in pacing_list]
# Starting at low frequency pacing (1Hz), 30 seconds at each pacing and then moving towards 230ms
for pacing in pacing_list:
    beats_per_pace = 30000/pacing
    p.schedule(1, start = offset, duration =0.5, period = pacing, multiplier = beats_per_pace)
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
start, duration, thresh = ap_duration(d, 30000*len(pacing_list), repolarisation = 95)

# First offset equal to zero, so remove first entry from the list
offset_list = offset_list[1:]

# Numpy array to contain final APD for each pacing cycle (and 3 previous APs in case of alternans)
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
    final_apd[i] = duration[index_start]
    final_apd2[i] = duration[index_start - 1]
    final_apd3[i] = duration[index_start - 2]
    final_apd4[i] = duration[index_start - 3]
    di.append(period[i]-duration[index_start -1])
    di2.append(period[i]-duration[index_start -2])
    di3.append(period[i]-duration[index_start -3])
    di4.append(period[i]-duration[index_start -4])

# The final peak doesn't have a peak after it, so can be indexed by the final duration recorded
final_apd[-1] = duration[-1]
di.append(period[-1] - duration[-2])

# If the graph has a long-short pattern, take previous APDs and DIs as well
final_apd2[-1] = duration[-2]
di2.append(period[-1] - duration[-3])
final_apd3[-1] = duration[-2]
di3.append(period[-1] - duration[-3])
final_apd4[-1] = duration[-2]
di4.append(period[-1] - duration[-3])


pl.figure()
pl.plot(d['engine.time'],d['membrane.V'])
pl.xlabel('Time (ms)')
pl.ylabel('')

# Plot the restitution curve
pl.figure()
pl.plot(pacing_list, final_apd, 'x', c = 'b')
pl.plot(pacing_list,final_apd2,'x', c = 'b')
pl.plot(pacing_list,final_apd3,'x', c = 'b')
pl.plot(pacing_list,final_apd4,'x', c = 'b')
pl.xlabel('PCL (ms)')
pl.ylabel('APD 95 (ms)')

pl.figure()
pl.plot(di, final_apd, 'x', c = 'b')
pl.plot(di2,final_apd2,'x', c = 'b')
pl.plot(di3,final_apd3,'x', c = 'b')
pl.plot(di4,final_apd4,'x', c = 'b')
pl.xlabel('DI (ms)')
pl.ylabel('APD 95 (ms)')
pl.xlim(0,360)
pl.ylim(190,270)

pl.title('Ohara (2011) Endothelial Cells Dynamic Protocol Restitution Curve')
pl.show()
