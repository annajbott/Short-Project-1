#!/usr/bin/env python

import matplotlib.pyplot as pl
import myokit
import numpy as np

def steady_state(m, pcl):
    p = myokit.pacing.blocktrain(pcl, 0.5, offset=0, level=1.0, limit=0)
    s = myokit.Simulation(m,p)
    '''
    d = s.run(pcl)
    value = max(d['membrane.V'])
    max_AP_time = d['membrane.V'].index(value)
    for i in range(max_AP_time + 1, len(d['membrane.V'])):
        if d['membrane.V'][i] <= value:
            value = d['membrane.V'][i]
        else:
            end_repo = i
            time_AP = round(d['engine.time'][i],0) + 1
            break

    #pl.figure()
    #pl.plot(d['engine.time'],d['membrane.V'])
    #pl.show()
    s.reset()
    print end_repo
    '''
    variation = np.zeros((40, 300))
    l = 0
    for i in range(10, 50):
        s.reset()
        d = s.run(i)
        variation[l][0] = (d['membrane.V'][-1])
        for j in range(1,20):
            d = s.run(pcl)
            variation[l][j] = (d['membrane.V'][-1])
        l += 1
    a = np.std(variation, axis = 1)
    time_AP = np.nonzero(a == max(a))[0][0] + 10
    print time_AP
    # Run simulation
    maxi = 600
    volt = np.zeros(maxi)

    s.reset()
    # One run at
    d = s.run(time_AP)
    volt[0] = (d['membrane.V'][-1])
    for i in range(1, maxi - 1):

      # Run a simulation of PCL length on top, so same point of AP is listed last
      d = s.run(pcl)

      # Append the list with the last recorded membrane potential
      volt[i] = (d['membrane.V'][-1])
      if np.all([volt[i-10:i] >= volt[i] - volt[i]*0.0001, volt[i-10:i] <= volt[i] + volt[i]*0.0001]) and i >= 10:
        steady = i
        break
    volt = volt[0:steady + 1]
    return(steady, volt)

def main():
    pcl = 1000
    m = myokit.load_model('grandi-2010.mmt')

    std, volt = steady_state(m, pcl)
    print std
    pl.figure()
    pl.plot(range(0,len(volt)), volt)
    pl.show()



if __name__ == "__main__":
    main()
