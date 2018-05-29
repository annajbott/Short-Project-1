#!/usr/bin/env python
import myokit
import matplotlib.pyplot as pl
from manual_APD import ap_duration

## Dynamic protocol diagram ##

m = myokit.load_model('ohara-2011.mmt')

# Set cell type
cell_types = {'Endocardial': 0, 'Epicardial' : 1, 'Mid-myocardial' : 2}
cell_type = 'Endocardial'

pacing_list = [600,300]
offset = 0
offset_list = []

p = myokit.Protocol()

for pacing in pacing_list:
    beats_per_pace = 30000/pacing
    if pacing == 300:
        beats_per_pace = 4
    p.schedule(1, start = offset, duration =0.5, period = pacing, multiplier = beats_per_pace)
    offset_list.append(offset)
    # Next set of pacing events to be scheduled by this offset (beats per pace * pace)
    offset += beats_per_pace*pacing


# Set up simulation using this scheduled protocol
s = myokit.Simulation(m, p)
s.set_constant('cell.mode', cell_types[cell_type])

# Run the simulation with final offset value, equal to time passed for whole protocol
d = s.run(offset, log = ['membrane.V','engine.time'])

# Use ap_duration function to calculate start times and durations
start, duration, thresh = ap_duration(d, 30000 + 3000, repolarisation = 90)

pl.figure()
pl.plot(d['engine.time'],d['membrane.V'])
pl.xlabel('Time (ms)')
pl.ylabel('Membrane Potential (mV)')

if len(start) > len(duration):
    start = start[0:-1]
for i, start in enumerate(start):
    duration_ap = duration[i]
    if i != 50 and i < 54 and i >47:
        pl.arrow(start, thresh[i], duration_ap, 0, head_width=4, head_length=80,
        length_includes_head=True, color = (0.1, 0.1, 0.1, 0.2))
        pl.text(start + 50, -72, str(int(duration_ap)) + ' ms')
    if i == 49:
        pacing_start = start
        pl.arrow(start, -95, 600, 0, head_width=5, head_length=100,
        length_includes_head=True, color = 'green')
        pl.text(start + 230, -100, 'PCL 1',  weight = 'bold', size = 'larger')


    if i == 50:
        pl.arrow(start, thresh[i], duration_ap, 0, head_width=4, head_length=80,
        length_includes_head=True, color = 'cyan')
        pl.text(start + 70, -67, 'APD',  weight = 'bold', size = 'larger')
        pl.text(start + 50, -72, str(int(duration_ap)) + ' ms')
        pl.arrow(start , -95, 300, 0, head_width=4, head_length=90,
        length_includes_head=True, color = 'magenta')
        pl.text(start + 70, -100, 'PCL 2',  weight = 'bold', size = 'larger')

pl.text(31240, -88, '...', weight = 'bold', size = 'x-large')
pl.xlim(28550, 31650)
pl.ylim(-105,43)

# Sub-plots
a = pl.axes([0.78, 0.6, 0.1, 0.15], facecolor='w')
pl.xticks([])
pl.yticks([])
y = [1000,1000,800,600,400,300,200,150,100,50]
x = [0,1,2,3,4,5,6,7,8,9]
pl.step(x,y)
pl.xlabel('Time')
pl.ylabel('PCL')
b = pl.axes([0.78, 0.4, 0.1, 0.15], facecolor='w')
pl.xticks([])
pl.yticks([])
y = [299.46968054, 297.91403479, 296.59887429, 295.32343115, 294.1510184, 292.849664, 291.51779892, 290.11774739, 288.45693228, 286.69384919, 284.78793966, 282.58996072, 280.23179626, 277.4542738, 274.42042784, 270.98010567, 267.28325386, 263.14119096, 258.51788763, 253.41075121, 247.98045541, 242.04006504, 235.49635441, 228.52240468, 221.05204964, 213.03323317, 204.56110355, 195.65527907, 186.49516962]
pl.plot([800, 780, 760, 740, 720, 700, 680, 660, 640, 620, 600, 580, 560, 540, 520, 500, 480, 460, 440, 420, 400, 380, 360, 340, 320, 300, 280, 260, 240],y)
pl.xlim(80, 820)
pl.ylim(150,330)
pl.xlabel('PCL')
pl.ylabel('APD')


pl.show()
