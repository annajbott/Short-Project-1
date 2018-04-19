#!/usr/bin/env python

import myokit
import matplotlib.pyplot as pl
from manual_APD import ap_duration
import numpy as np

# Get the model
m = myokit.load_model('tentusscher-2006.mmt')

# Initial pacing, step size to decrease and minimum to decrease to
pcl = 800
step_size = -20
min_pcl = 40
paces = (pcl - min_pcl)/abs(step_size) + 1

# How many beats to do per pacing cycle length
beats_per_pace = 50

# Pre pace at this intial pacing
p = myokit.pacing.blocktrain(pcl, 0.5, offset=0, level=1.0, limit=0)
s = myokit.Simulation(m, p,apd_var='membrane.V')
s.pre(pcl * 40)

# Set threshold for AP repolarisation
threshold = 0.9 * s.state()[m.get('membrane.V').indice()]

# Empty protocol
p = myokit.Protocol()
s = myokit.Simulation(m, p)
s.reset()

# Empty arrays to fill.
period = []
offset_list = []

# Start with offset in schedule equal to 0, offset delays when to schedule next PCL change
offset = 0

# Decrease PCL until it hits the minimum value, scheduling new pacing events
while pcl > min_pcl:

  #Initialising protocol with new pacing
  p.schedule(level = 1, start = offset, duration = 0.5, period = pcl, multiplier = beats_per_pace )
  offset_list.append(offset)

  # Next set of pacing events to be scheduled by this offset (beats per pace * pace)
  offset += beats_per_pace*pcl
  period.append(pcl)
  pcl += step_size

# Set up simulation using this scheduled protocol
s = myokit.Simulation(m, p)

# Run the simulation with final offset value, equal to time passed for whole protocol
d = s.run(offset)

# Use ap_duration function to calculate start times and durations
start, duration = ap_duration(d, paces*beats_per_pace,threshold)

# First offset equal to zero, so remove first entry from the list
offset_list = offset_list[1:]

# Numpy array to contain final APD for each pacing cycle
final_apd = np.zeros(len(offset_list) + 1)

print offset_list
print start

# Fill final_apd array
for i in range(len(offset_list)):
    # Final peak of pacing cycle = peak before the first of a new pacing cycle
    # Index_start = The index of the start of pacing cycle in start array

    index_start = np.nonzero(np.round(start,0) >= offset_list[i])[0][0]

    # index_1 -1 to get peak at the end of the previous cycle
    final_apd[i] = duration[index_start - 1]

# The final peak doesn't have a peak after it, so can be indexed by the final duration recorded
final_apd[-1] = duration[-1]

# Plot the restitution curve
pl.figure()
pl.plot(period, final_apd, 'x')
pl.show()
