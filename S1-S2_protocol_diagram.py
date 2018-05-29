#!/usr/bin/env python
import myokit
import matplotlib.pyplot as pl
from manual_APD import ap_duration

## S1-S2 protocol diagram ##

m = myokit.load_model('tentusscher-2006.mmt')

# Pacing for S1 beats is 600ms
pcl = 600
p = myokit.pacing.blocktrain(pcl, 0.5, offset=0, level=1.0, limit=0)
s = myokit.Simulation(m, p)
s.pre(pcl*50)
s.set_constant('cell.type', 1)
p = myokit.Protocol()
p.schedule(1,100,0.5,pcl,3)
p.schedule(1,2*pcl + 500, 0.5, pcl, 1)
s.set_protocol(p)
d = s.run(4*pcl)
start, duration, thresh = ap_duration(d, 4)

pl.figure()
pl.plot(d['engine.time'],d['membrane.V'])
if len(start) > len(duration):
    start = start[0:-1]
for i, start in enumerate(start):
    duration_ap = duration[i]
    if i< 3:
        pl.arrow(start, thresh[i], duration_ap, 0, head_width=4, head_length=80,
        length_includes_head=True, color = (0.1, 0.1, 0.1, 0.2))
        pl.text(start + 50, -85, str(int(duration_ap)) + ' ms')
    if i == 2:
        DI_start = start + duration_ap
        tr = thresh[i]
    if i == 3:
        pl.arrow(start, thresh[i], duration_ap, 0, head_width=5, head_length=100,
        length_includes_head=True, color = 'cyan')
        pl.text(start + 50, -70, 'APD',  weight = 'bold', size = 'larger')
        pl.text(start + 50, -85, str(int(duration_ap)) + ' ms')
        pl.text(start + 30, 23, 'S2 beat', weight = 'bold', size = 'x-large')
        pl.arrow(DI_start, tr, start - DI_start, 0,head_width=3, head_length=50,
          length_includes_head=True, color = 'red')
        pl.text(DI_start + 25, -70, 'DI',  weight = 'bold', size = 'larger')
        pl.text(DI_start + 16, -80, str(int(start - DI_start)) + ' ms')


pl.xlim(0,2400)
pl.ylim(-91,40)
pl.xlabel('Time (ms)')
pl.ylabel('Membrane Potential (mV)')
a = pl.axes([0.78, 0.6, 0.1, 0.15], facecolor='w')
pl.plot([66.338310742602289, 91.333092543839882, 116.34367249695111, 141.33883632529449, 166.31995200959534, 191.32480735791501, 216.32133758529778, 241.32932267397274, 266.32631038000437, 291.33831074260229, 316.33091475695875, 341.34364864482922, 366.32775930255212, 391.33338327761294, 416.32338107664538, 441.32956933791138, 466.32868093508205, 491.31388944065185, 516.32932267397268, 541.33090774378911, 566.36249883712492, 591.32452892858419, 616.3326298874897, 641.32193863542011, 666.32429270043542, 691.33710974555493],[207.21279404191733, 226.24513291490268, 241.31621082692033, 253.16235805577173, 262.42073561522079, 269.73559263607234, 275.49360778728021, 280.09376979593787, 283.80112175048907, 286.79272617334914, 289.18920380236398, 291.20297944234449, 292.88901832174906, 294.20379461035759, 295.35913852281851, 296.39503742508305, 297.14626974769089, 297.84380598246491, 298.41137025941021, 298.90532929562994, 299.33394047169782, 299.69936657175208, 299.97529063953277, 300.24654019795616, 300.41777356035084, 300.65193021772586])
pl.xticks([])
pl.yticks([])
pl.xlabel('DI (ms)')
pl.ylabel('APD 90 (ms)')
pl.show()
