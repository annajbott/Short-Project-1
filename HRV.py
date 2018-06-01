#!/usr/bin/env python
import myokit
import matplotlib.pyplot as pl
from manual_APD import ap_duration
import numpy as np
from HF_model import Ord_HF_Gomez
from restitution_protocols import *


# Protocol
p = myokit.Protocol()
hf_p = myokit.Protocol()
hf_6minwalk_p = myokit.Protocol()
m = myokit.load_model('ohara-2011.mmt')
hf_m = Ord_HF_Gomez()
number_runs = 50
# HRV
# Need range of HR -- amplitude (like 50bpm to 95bpm)
# Need frequency, slower breathing--> smaller frequency (longer period)

# 0.15-0.4Hz breathing rate (9-24 breaths a minute / 2500ms - 6666ms)
# Make protocol with HRV included first for no HF (50-90 pattern)

# Using sin wave to simulate step protocol
points = np.linspace(0,2*np.pi, num  = 30, endpoint = False)
hr_pacing_list = []
hf_hr_pacing_list = []
hf_6minwalk_pacing_list = []
average = 70
hf_average = 90
# Weighted average from Correa, Silva Alves, Bianchim (2013)
hf_6minwalk_average = 114
for i in points:
    # Using a sin wave centered around the average (18.5,10, 5 = amplitude)
    hr_pacing_list.append(average +18.5*np.sin(i-np.pi/2.0))
    hf_hr_pacing_list.append(hf_average +10*np.sin(i-np.pi/2.0))
    hf_6minwalk_pacing_list.append(hf_6minwalk_average +8*np.sin(i-np.pi/2.0))

# Numpy arrays
hr_pacing_list = np.asarray(hr_pacing_list)
hf_hr_pacing_list = np.asarray(hf_hr_pacing_list)
hf_6minwalk_pacing_list = np.asarray(hf_6minwalk_pacing_list)

# Using guestimate from points on graphic (non HF)
#hr_pacing_list = np.array([52,52.5,55,60,65,74,80,82.5,85,88,89,87,70.5])
#hf_hr_pacing_list = np.array([80, 80.27, 81.62, 84.32,87.03, 91.89, 95.13, 96.48, 97.83, 99.45, 99.99, 98.91,89.99])

# Generate integer pacing values from instantaneous HR
pacing_list = [np.round(60*1000.0/element,0) for element in hr_pacing_list]
hf_pacing_list = [np.round(60*1000.0/element,0)  for element in hf_hr_pacing_list]
exercise_pacing_list = [np.round(60*1000.0/element,0)  for element in hf_6minwalk_pacing_list]
beats_per_pace = 1

# Setting initital offset to zero
offset = 0
hf_offset = 0
hf_6minwalk_offset = 0

# Empty array to store offsets of new pacing in
offset_array = []
hf_offset_array = []
hf_6minwalk_offset_array =[]

# Array to store every PCL in
pace_array = []
hf_pace_array = []
hf_6minwalk_pace_array = []

for number in range(0, number_runs):
    # Create protocol for healthy heart pacing
    for pacing in pacing_list:
        p.schedule(1, start = offset, duration = 0.5, period = pacing, multiplier = beats_per_pace)
        offset += beats_per_pace*pacing
        pace_array.append(pacing)
        offset_array.append(offset)
    # Create protocol for HF at rest
    for pacing in hf_pacing_list:
        hf_p.schedule(1, start = hf_offset, duration = 0.5, period = pacing, multiplier = beats_per_pace)
        hf_offset += beats_per_pace*pacing
        hf_pace_array.append(pacing)
        hf_offset_array.append(hf_offset)
    # Create protocol for HF during final minutes of 6-min walk test
    for pacing in exercise_pacing_list:
        hf_6minwalk_p.schedule(1, start = hf_6minwalk_offset, duration = 0.5, period = pacing, multiplier = beats_per_pace)
        hf_6minwalk_offset += beats_per_pace*pacing
        hf_6minwalk_pace_array.append(pacing)
        hf_6minwalk_offset_array.append(hf_6minwalk_offset)


## Create and run simulation for healthy ventricular cell model ##
## ------------------------------------------------------------ ##

s = myokit.Simulation(m, p)
s.pre(np.sum(pacing_list)*10)
s.reset()
d = s.run(offset, log = ['membrane.V','engine.time'])

# Plot APs for healthy ventricular cell model
pl.figure()
pl.plot(d['engine.time'],d['membrane.V'])
pl.xlabel('Time (ms)')
pl.ylabel('Membrane Potential (mV)')
# Calculate APDs
start, duration, thresh = ap_duration(d, number_runs*len(pacing_list), repolarisation = 90)

for i in range(0, 21):
    # Plot arrows and PCL for first sin wave of APs
    pl.arrow(start[i], thresh[i], start[i+1]-start[i], 0, head_width=3, head_length=80,
    length_includes_head=True, color = 'green')
    pl.text(start[i] + 20, -72, str(pace_array[i]))

## Create and run simulation for HF at rest ##
## ---------------------------------------- ##

hf_s = myokit.Simulation(hf_m, hf_p)
hf_s.pre(np.sum(hf_pacing_list)*10)
hf_s.reset()
hf_d = hf_s.run(hf_offset, log = ['membrane.V','engine.time'])

# Plot APs for HF ventricular cell model at rest
pl.figure()
pl.plot(hf_d['engine.time'],hf_d['membrane.V'])
pl.title('HF')
# Calculate APDs
hf_start, hf_duration, hf_thresh = ap_duration(hf_d, number_runs*len(pacing_list), repolarisation = 90)

for i in range(0, 21):
    # Plot arrows and PCL for first sin wave of APs
    pl.arrow(hf_start[i], hf_thresh[i], hf_start[i+1]-hf_start[i], 0, head_width=3, head_length=80,
    length_includes_head=True, color = 'green')
    pl.text(hf_start[i] + 20, -72, str(hf_pace_array[i]))


## Create and run simulation for HF during final minutes of 6 min walk test ##
## ------------------------------------------------------------------------ ##

hf_6minwalk_s = myokit.Simulation(hf_m, hf_6minwalk_p)
hf_6minwalk_s.pre(np.sum(exercise_pacing_list)*10)
hf_6minwalk_s.reset()
hf_6minwalk_d = hf_6minwalk_s.run(hf_6minwalk_offset, log = ['membrane.V','engine.time'])

# Plot APs for HF ventricular cell model 6 min walk
pl.figure()
pl.plot(hf_6minwalk_d['engine.time'],hf_6minwalk_d['membrane.V'])
pl.title('HF 6-min walk test')
# Calculate APDs
hf_6minwalk_start, hf_6minwalk_duration, hf_6minwalk_thresh = ap_duration(hf_6minwalk_d, number_runs*len(pacing_list), repolarisation = 90)


for i in range(0, 21):
    # Plot arrows and PCL for first sin wave of APs
    pl.arrow(hf_6minwalk_start[i], hf_6minwalk_thresh[i], hf_6minwalk_start[i+1]-hf_6minwalk_start[i], 0, head_width=3, head_length=80,
    length_includes_head=True, color = 'green')
    pl.text(hf_6minwalk_start[i] + 20, -72, str(hf_6minwalk_pace_array[i]))
    if  i == 20:
        pl.arrow(hf_6minwalk_start[i], hf_6minwalk_thresh[i] - 15, hf_6minwalk_duration[i], 0, head_width=3, head_length=80,
        length_includes_head=True, color = 'red')
        pl.text(hf_6minwalk_start[i] + 20, -80, str(hf_6minwalk_duration[i]))




# Figure to plot HR (BPM) against time
hr_bpm = [60*1000.0/element for element in pace_array]
hf_hr_bpm = [60*1000.0/element for element in hf_pace_array]
hf_6minwalk_bpm = [60*1000.0/element for element in hf_6minwalk_pace_array]


pl.figure()
pl.plot(offset_array,hr_bpm,'.-')
pl.plot(hf_offset_array, hf_hr_bpm,'.-')
pl.plot(hf_6minwalk_offset_array, hf_6minwalk_bpm,'.-')
pl.xlim(0,50000)
pl.xlabel('Time (ms)')
pl.ylabel('Instantaneous HR (BPM)')
pl.title('Dynamic Pacing Protocol to mimic HRV')
pl.legend(['Healthy HRV', 'HF HRV at rest', 'HF HRV 6-min walk test'])


pl.figure()
#pace_array = np.repeat(pace_array,beats_per_pace)
pl.plot(pace_array[len(pace_array)/2:], duration[len(pace_array)/2:], '.')
pl.plot(hf_pace_array[len(pace_array)/2:], hf_duration[len(pace_array)/2:], '.')
pl.plot(hf_6minwalk_pace_array[len(pace_array)/2 -1:-1], hf_6minwalk_duration[len(pace_array)/2:], '.')
pl.xlabel('PCL (ms)')
pl.ylabel('APD (ms)')
pl.title('Restitution Curves for Dynamic Protocol based on HRV')
pl.legend(['Healthy', 'HF at rest', 'HF 6-min walk test'])
pl.show()
