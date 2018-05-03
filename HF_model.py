#!/usr/bin/env python

import matplotlib.pyplot as pl
import myokit
import numpy as np
from manual_APD import ap_duration

### HF model O'hara 2011 ###

def Ord_HF_Gomez():
    m = myokit.load_model('ohara-2011.mmt')

    # Late sodium current- ORd model 180% in HF
    inal = m.get('inal.INaL')
    inal.set_rhs(myokit.Multiply(myokit.Number(1.8), inal.rhs()))

    # Time constant of inactivation of the INaL- ORd 180%
    # phosphorylated or non-phosphorylated?? th vs thp or both
    th = m.get('inal.th')
    th.set_rhs(myokit.Multiply(myokit.Number(1.8), th.rhs()))
    #thp = m.get('inal.thp')
    #thp.set_rhs(myokit.Multiply(myokit.Number(1.8), thp.rhs()))

    # Transient outward K+ current- ORd 40% in HF
    ito = m.get('ito.Ito')
    ito.set_rhs(myokit.Multiply(myokit.Number(0.4), ito.rhs()))

    # Inward rectifier K+ current- ORd 68% in HF
    ik1 = m.get('ik1.IK1')
    ik1.set_rhs(myokit.Multiply(myokit.Number(0.68), ik1.rhs()))

    # Na+/ K+ pump current- ORd 70% in HF
    inak = m.get('inak.INaK')
    inak.set_rhs(myokit.Multiply(myokit.Number(0.7), inak.rhs()))

    # Sodium Calcium Exchanger- ORd 175% in HF
    incx = m.get('inaca.INaCa')
    incx.set_rhs(myokit.Multiply(myokit.Number(1.75), incx.rhs()))

    # SERCA pump - ORd 50% in HF
    serc = m.get('serca.Jup') # Could be . Jtr or others, not clear
    serc.set_rhs(myokit.Multiply(myokit.Number(0.5), serc.rhs()))

    # Ileak - ORd 130% in HF
    leak = m.get('serca.Jleak')
    leak.set_rhs(myokit.Multiply(myokit.Number(1.3), leak.rhs()))

    # CamKa - ORd 150% in HF
    camka = m.get('camk.CaMKa')
    camka.set_rhs(myokit.Multiply(myokit.Number(1.5), camka.rhs()))

    # Jrel- Non-phosphorylated Ca2+ release via Ryr - ORd 80% in HF
    ryr = m.get('ryr.Jrel') # Think leak occurs via Ryanodine receptors, could be serca.Jleak
    ryr.set_rhs(myokit.Multiply(myokit.Number(0.8), ryr.rhs()))
    return(m)

def GPB_HF_Gomez():
    m = myokit.load_model('grandi-2010.mmt')

    # Late Na+ current - fast component where is late in model
    inal = m.get('inal.GNal') # GNal multiplied by both junc and sl
    inal.set_rhs(myokit.Multiply(myokit.Number(2.0), inal.rhs()))

    # Transient outward K+ current- GPB 40% in HF
    ito = m.get('ito.ito')
    ito.set_rhs(myokit.Multiply(myokit.Number(0.4), ito.rhs()))

    # Inward rectifier K+ current- GPB 68% in HF
    ik1 = m.get('ik1.I_k1')
    ik1.set_rhs(myokit.Multiply(myokit.Number(0.68), ik1.rhs()))

    # Na+/ K+ pump current- GPB 90% in HF
    inak = m.get('inak.I_nak')
    inak.set_rhs(myokit.Multiply(myokit.Number(0.9), inak.rhs()))

    # Background Na+ current- GPB 0% in HF
    inab = m.get('inab.GNaB') # Then slow and junc = 0, as both multiplied by GNaB
    inab.set_rhs(myokit.Multiply(myokit.Number(0.0), inab.rhs()))

    # Background Ca2+ current- GPB 153% in HF
    icab = m.get('icabk.GCaB') # Slow and junc, both multiplied by GCaB (ask anna/michael)
    icab.set_rhs(myokit.Multiply(myokit.Number(1.53), icab.rhs()))

    # Sodium/calcium exchanger current- GPB 175% in HF
    incx = m.get('incx.IbarNCX') # Both slow and junc multiplied
    incx.set_rhs(myokit.Multiply(myokit.Number(1.75), incx.rhs()))

    # Sarco/ endoplasmic reticulum Ca2+ pump current- GPB 50% in HF
    serca = m.get('caflux.J_serca')
    serca.set_rhs(myokit.Multiply(myokit.Number(0.5), serca.rhs()))

    # SR Ca2+ leak GPB 300% in HF
    leak = m.get('caflux.J_SRleak')
    leak.set_rhs(myokit.Multiply(myokit.Number(3.0), leak.rhs()))

    # SR Ca2+ leak GPB 300% in HF
    ec50 = m.get('caflux.ec50SR')
    ec50.set_rhs(myokit.Multiply(myokit.Number(0.89), ec50.rhs()))

    return(m)


def main():
    ## Basic cycle length and number of paces to run for
    bcl = 1000
    paces = 1

    # Set protocol
    p = myokit.pacing.blocktrain(bcl, 0.5, offset=100, level=1.0, limit=0)

    ## Grandi model ##
    ## ------------ ##

    # Set cell type
    cell_types = {'Endocardial' : 0, 'Epicardial' : 1} #'Mid-myocardial' : 2}

    pl.figure()
    for i, cell_type in enumerate(cell_types):
        m = myokit.load_model('grandi-2010.mmt')
        s = myokit.Simulation(m,p)

        s.set_constant('type.epi', cell_types[cell_type])
        s.pre(200*bcl)
        d = s.run(paces*bcl)
        pl.subplot(1, 2, cell_types[cell_type] + 1)
        pl.plot(d['engine.time'],d['membrane.V'])

        # Use Function to set parameters
        m = GPB_HF_Gomez()
        s = myokit.Simulation(m,p)

        s.set_constant('type.epi', cell_types[cell_type])
        s.pre(200*bcl)
        d1 = s.run(paces*bcl)
        pl.plot(d1['engine.time'],d1['membrane.V'])
        pl.legend(['Normal','HF '] )
        pl.ylabel('Membrane Potential (mV)')
        pl.xlabel('Time (ms)')
        pl.title("{} cell".format(cell_type))
    pl.suptitle("Grandi (2010) model at {} ms pace cycle length".format(bcl))

    ## O'hara model ##
    ## ------------ ##

    # Set cell type
    cell_types = {'Endocardial' : 0, 'Epicardial' : 1, 'Mid-myocardial' : 2}

    pl.figure()
    for i, cell_type in enumerate(cell_types):
        m = myokit.load_model('ohara-2011.mmt')
        s = myokit.Simulation(m,p)

        s.set_constant('cell.mode', cell_types[cell_type])
        s.pre(200*bcl)
        d = s.run(paces*bcl)
        start, duration, thresh = ap_duration(d, paces)
        pl.subplot(1, 3, cell_types[cell_type] + 1)
        pl.plot(d['engine.time'],d['membrane.V'])
        pl.text(500,0,"Healthy APD- {} ms".format(np.round(duration[0],2)))


        # Use Function to set parameters
        m = Ord_HF_Gomez()
        s = myokit.Simulation(m,p)

        s.set_constant('cell.mode', cell_types[cell_type])
        s.pre(200*bcl)
        d1 = s.run(paces*bcl)
        start, duration_hf, thresh = ap_duration(d1, paces)
        pl.plot(d1['engine.time'],d1['membrane.V'])
        pl.text(500,-5,"HF APD- {} ms".format(np.round(duration_hf[0],2)))
        pl.legend(['Normal','HF '] )
        pl.ylabel('Membrane Potential (mV)')
        pl.xlabel('Time (ms)')
        pl.title("{} cell".format(cell_type))
    pl.ylim(-93, 48)
    pl.suptitle("O'hara (2011) model at {} ms pace cycle length".format(bcl))
    pl.show()

if __name__ == "__main__":
    main()
