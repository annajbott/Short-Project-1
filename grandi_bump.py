#!/usr/bin/env python
import myokit
import matplotlib.pyplot as pl
from manual_APD import ap_duration
import numpy as np

## Investigating Grandi bump ##
## ------------------------- ##

m = myokit.load_model('grandi-2010.mmt')
cell_types = {'Endocardial': 0, 'Epicardial' : 1}
cell_type = 'Endocardial'

pcl = 1000
p = myokit.pacing.blocktrain(pcl, 0.5, offset=20, level=1.0, limit=0)
s = myokit.Simulation(m,p)
s.set_constant('type.epi', cell_types[cell_type])
s.set_tolerance(1*np.e**-8, 1*np.e**-8)

d = s.run(pcl)
pl.figure()
pl.subplot(3,2,1)
pl.plot(d['engine.time'],d['membrane.V'])
pl.title('AP')
pl.subplot(3,2,2)
pl.plot(d['engine.time'],d['ito.ito'])
pl.title('Ito')
pl.subplot(3,2,3)
pl.plot(d['engine.time'],d['ikr.I_kr'])
pl.title('IKr')
pl.subplot(3,2,4)
pl.plot(d['engine.time'],d['ik1.I_k1'])
pl.title('IK1')
pl.subplot(3,2,5)
pl.plot(d['engine.time'],d['ical.I_Ca_junc'])
pl.plot(d['engine.time'],d['ical.I_Ca_sl'])
pl.legend(['Junctional cleft facing', 'Subsarcolemmal space injecting'])
pl.title('ICal')
pl.subplot(3,2,6)
pl.plot(d['engine.time'],d['calcium.I_Ca_tot_junc'])
pl.plot(d['engine.time'],d['calcium.I_Ca_tot_sl'])
pl.legend(['Junctional cleft facing', 'Subsarcolemmal space injecting'])
pl.title('I Ca')
pl.suptitle('{} 1st beat, Grandi (2010)'.format(cell_type))

# Reset and run 200 beats at 1000ms
s.reset()
s.set_constant('type.epi', cell_types[cell_type])
s.set_tolerance(1*np.e**-8, 1*np.e**-8)
s.pre(pcl*199)
d = s.run(pcl)
pl.figure()
pl.subplot(3,2,1)
pl.plot(d['engine.time'],d['membrane.V'])
pl.title('AP')
pl.subplot(3,2,2)
pl.plot(d['engine.time'],d['ito.ito'])
pl.title('Ito')
pl.subplot(3,2,3)
pl.plot(d['engine.time'],d['ikr.I_kr'])
pl.title('IKr')
pl.subplot(3,2,4)
pl.plot(d['engine.time'],d['ik1.I_k1'])
pl.title('IK1')
pl.subplot(3,2,5)
pl.plot(d['engine.time'],d['ical.I_Ca_junc'])
pl.plot(d['engine.time'],d['ical.I_Ca_sl'])
pl.legend(['Junctional cleft facing', 'Subsarcolemmal space injecting'])
pl.title('ICal')
pl.subplot(3,2,6)
pl.plot(d['engine.time'],d['calcium.I_Ca_tot_junc'])
pl.plot(d['engine.time'],d['calcium.I_Ca_tot_sl'])
pl.legend(['Junctional cleft facing', 'Subsarcolemmal space injecting'])
pl.title('I Ca')
pl.suptitle('{} 200th beat Grandi (2010)'.format(cell_type))


# Reset and run 200 beats at 1000ms
s = myokit.Simulation(m,p)
s.set_constant('type.epi', cell_types[cell_type])
s.set_tolerance(1*np.e**-8, 1*np.e**-8)

pl.figure()
beats = [1,2,3,6,12,20,60,100,150,200]
time_beats_start = [(element -1)*10000 for element in beats]
beats_legend = [str(element)+'th AP' for element in beats]
beats_legend[0] = '1st AP'
beats_legend[1] = '2nd AP'
beats_legend[2] = '3rd AP'

#num_colours = len(beats)
#cm = pl.get_cmap('gist_rainbow')
#pl.set_color_cycle([cm(1.*i/num_colours) for i in range(num_colours)])

d = s.run(200*pcl,log_interval = 0.1, log = ['engine.time','membrane.V','calcium.I_Ca_tot_junc','calcium.I_Ca_tot_sl', 'ical.I_Ca_junc', 'ical.I_Ca_sl' ])
pl.subplot(1,5,1)
for beat in time_beats_start:
    pl.plot(d['engine.time'][0:10000],d['membrane.V'][beat:beat+10000])
pl.title('Membrane potential, AP')
pl.subplot(1,5,2)
for beat in time_beats_start:
    pl.plot(d['engine.time'][0:10000],d['ical.I_Ca_junc'][beat:beat+10000])
pl.title('ICal (junc)')
pl.subplot(1,5,3)
for beat in time_beats_start:
    pl.plot(d['engine.time'][0:10000],d['ical.I_Ca_sl'][beat:beat+10000])
pl.title('ICal (sl)')
pl.subplot(1,5,4)
for beat in time_beats_start:
    pl.plot(d['engine.time'][0:10000],d['calcium.I_Ca_tot_junc'][beat:beat+10000])
pl.title('I Ca (junc)')
pl.subplot(1,5,5)
for beat in time_beats_start:
    pl.plot(d['engine.time'][0:10000],d['calcium.I_Ca_tot_sl'][beat:beat+10000])
    pl.legend(beats_legend)
pl.title('I Ca (sl)')
pl.suptitle('Membrane potential and calcium channel currents at a certain beat in a cycle using the Grandi (2010) Model')






## Checking calcium channels at specific point across 500 APs
s.reset()
s = myokit.Simulation(m,p)
s.set_constant('type.epi', cell_types[cell_type])
s.set_tolerance(1*np.e**-8, 1*np.e**-8)
d = s.run(500*pcl, log_interval = 0.1, log = ['engine.time','calcium.I_Ca_tot_junc','calcium.I_Ca_tot_sl', 'ical.I_Ca_junc', 'ical.I_Ca_sl' ])
pl.figure()
point_time = d['engine.time'][650::10000]
point_AP = [element/1000 for element in point_time]
cajunc_time = d['calcium.I_Ca_tot_junc'][650::10000]
casl_time = d['calcium.I_Ca_tot_sl'][650::10000]
icaljunc_time = d['ical.I_Ca_junc'][650::10000]
icalsl_time = d['ical.I_Ca_sl'][650::10000]

pl.subplot(2,2,1)
pl.plot(point_AP,cajunc_time)
pl.xlabel('AP')
pl.ylabel('Channel Current (mV)')
pl.title('I Ca junc')

pl.subplot(2,2,2)
pl.plot(point_AP, casl_time)
pl.xlabel('Time (ms)')
pl.ylabel('Channel Current (mV)')
pl.title('I Ca Sl')
pl.suptitle('Ca current at 65ms into AP')

pl.subplot(2,2,3)
pl.plot(point_AP, icaljunc_time)
pl.xlabel('AP')
pl.ylabel('Channel Current (mV)')
pl.title('ICal junc')
pl.suptitle('Ca current at 65ms into AP')

pl.subplot(2,2,4)
pl.plot(point_AP, icalsl_time)
pl.xlabel('AP')
pl.ylabel('Channel Current (mV)')
pl.title('ICal Sl')
pl.suptitle('Ca current at 65ms into AP, with high precision solver')

pl.show()
