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


s.reset()
s.set_constant('type.epi', cell_types[cell_type])
d = s.run(500*pcl, log_interval = 0.1, log = ['engine.time','calcium.I_Ca_tot_junc','calcium.I_Ca_tot_sl', 'ical.I_Ca_junc', 'ical.I_Ca_sl' ])
pl.figure()
point_time = d['engine.time'][450::10000]

cajunc_time = d['calcium.I_Ca_tot_junc'][450::10000]
casl_time = d['calcium.I_Ca_tot_sl'][450::10000]
icaljunc_time = d['ical.I_Ca_junc'][450::10000]
icalsl_time = d['ical.I_Ca_sl'][450::10000]

pl.subplot(2,2,1)
pl.plot(point_time,cajunc_time)
pl.xlabel('Time (ms)')
pl.ylabel('Membrane Potential (mV)')
pl.title('I Ca junc')

pl.subplot(2,2,2)
pl.plot(point_time, casl_time)
pl.xlabel('Time (ms)')
pl.ylabel('Membrane Potential (mV)')
pl.title('I Ca Sl')
pl.suptitle('Ca current at 45ms into AP')

pl.subplot(2,2,3)
pl.plot(point_time, icaljunc_time)
pl.xlabel('Time (ms)')
pl.ylabel('Membrane Potential (mV)')
pl.title('ICal junc')
pl.suptitle('Ca current at 45ms into AP')

pl.subplot(2,2,4)
pl.plot(point_time, icalsl_time)
pl.xlabel('Time (ms)')
pl.ylabel('Membrane Potential (mV)')
pl.title('ICal Sl')
pl.suptitle('Ca current at 45ms into AP')

pl.show()
