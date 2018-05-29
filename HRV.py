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
m = myokit.load_model('ohara-2011.mmt')
hf_m = Ord_HF_Gomez()
# HRV
# Need range of HR -- amplitude (like 50bpm to 95bpm)
# Need frequency, slower breathing--> smaller frequency (longer period)

# 0.15-0.4Hz breathing rate (9-24 breaths a minute / 2500ms - 6666ms)
# Make protocol with HRV included first for no HF (50-90 pattern)


# Using guestimate from points on graphic (non HF)
hr_pacing_list = np.array([52,52.5,55,60,65,74,80,82.5,85,88,89,87,70.5])
hf_hr_pacing_list = np.array([80, 80.27, 81.62, 84.32,87.03, 91.89, 95.13, 96.48, 97.83, 99.45, 99.99, 98.91,89.99])
pacing_list = [int(60*1000/element) for element in hr_pacing_list]
hf_pacing_list = [int(60*1000/element) for element in hf_hr_pacing_list]
beats_per_pace = 1
offset = 0
offset_array = []
pace_array = []

hf_offset = 0
hf_offset_array = []
hf_pace_array = []
for number in range(0, 100):
    for pacing in pacing_list:
        p.schedule(1, start = offset, duration = 0.5, period = pacing, multiplier = beats_per_pace)
        offset += beats_per_pace*pacing
        pace_array.append(pacing)
        offset_array.append(offset)
    for pacing in hf_pacing_list:
        hf_p.schedule(1, start = hf_offset, duration = 0.5, period = pacing, multiplier = beats_per_pace)
        hf_offset += beats_per_pace*pacing
        hf_pace_array.append(pacing)
        hf_offset_array.append(hf_offset)




s = myokit.Simulation(m, p)
s.pre(np.sum(pacing_list)*10)
s.reset()
d = s.run(offset, log = ['membrane.V','engine.time'])

s1 = myokit.Simulation(hf_m, hf_p)
s1.pre(np.sum(pacing_list)*10)
s1.reset()
d1 = s1.run(offset, log = ['membrane.V','engine.time'])


pl.figure()
pl.plot(d['engine.time'],d['membrane.V'])
pl.xlabel('Time (ms)')
pl.ylabel('Membrane Potential (mV)')
pl.figure()
pl.plot(d1['engine.time'],d1['membrane.V'])
pl.title('HF')

start, duration, thresh = ap_duration(d, 100*len(pacing_list), repolarisation = 90)
hf_start, hf_duration, hf_thresh = ap_duration(d1, 100*len(pacing_list), repolarisation = 90)


'''
for i in range(90, 110):

    pl.arrow(start[i], -90, start[i+1]-start[i], 0, head_width=3, head_length=80,
    length_includes_head=True, color = 'green')
    pl.text(start[i] + 20, -72, str(int(start[i+1]-start[i])) + ' pacing = {}'.format(pace_array[i]))

'''
# Figure to plot HR (BPM) against
hr_bpm = [60*1000.0/element for element in pace_array]
hf_hr_bpm = [60*1000.0/element for element in hf_pace_array]

pl.figure()
pl.plot(offset_array,hr_bpm)
pl.plot(hf_offset_array, hf_hr_bpm)

pl.figure()
#pace_array = np.repeat(pace_array,beats_per_pace)
pl.plot(pace_array[len(pace_array)/2:], duration[len(pace_array)/2:], 'x')
pl.plot(hf_pace_array[len(pace_array)/2:], hf_duration[len(pace_array)/2:], 'o')

pl.show()
