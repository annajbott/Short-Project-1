#!/usr/bin/env python

import matplotlib.pyplot as pl
import myokit
import numpy as np
from manual_APD import ap_duration

### HF model O'hara 2011 ###

def Ord_HF_Gomez(cell_type = None):
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
    # If no cell type passed, +75% for all cell types
    if cell_type == None:
        incx.set_rhs(myokit.Multiply(myokit.Number(1.75), incx.rhs()))
    # Specified to be epicardial cell (1)
    elif cell_type == 1:
        incx.set_rhs(myokit.Multiply(myokit.Number(2.0), incx.rhs()))
    # Endocardial or mid-myocardial (0 or 2)
    else:
        incx.set_rhs(myokit.Multiply(myokit.Number(1.6), incx.rhs()))

    # SERCA pump - ORd 50% in HF
    serc = m.get('serca.Jup') # Could be . Jtr or others, not clear
    # No cell type specified, set all cell types the same (-50%)
    if cell_type == None:
        serc.set_rhs(myokit.Multiply(myokit.Number(0.5), serc.rhs()))
    # Epicardial
    elif cell_type == 1:
        serc.set_rhs(myokit.Multiply(myokit.Number(0.75), serc.rhs()))
    # Mid-myocardial
    elif cell_type == 2:
        serc.set_rhs(myokit.Multiply(myokit.Number(0.6), serc.rhs()))
    # Endocardial
    else:
        serc.set_rhs(myokit.Multiply(myokit.Number(0.45), serc.rhs()))

    # Ileak - ORd 130% in HF
    leak = m.get('serca.Jleak')
    leak.set_rhs(myokit.Multiply(myokit.Number(1.3), leak.rhs()))

    # CamKa - ORd 150% in HF
    camka = m.get('camk.CaMKa')
    camka.set_rhs(myokit.Multiply(myokit.Number(1.5), camka.rhs()))

    # Jrel- Non-phosphorylated Ca2+ release via Ryr - ORd 80% in HF
    ryr = m.get('ryr.Jrel')
    ryr.set_rhs(myokit.Multiply(myokit.Number(0.8), ryr.rhs()))
    return(m)

def GPB_HF_Gomez(cell_type = None):
    m = myokit.load_model('grandi-2010.mmt')

    # Original grandi model does not contain late channel. Changes from Tenor et al.
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
    if cell_type == None:
        incx.set_rhs(myokit.Multiply(myokit.Number(1.75), incx.rhs()))
    # Grandi model no mid-myocardial cells
    elif cell_type == 0:
        incx.set_rhs(myokit.Multiply(myokit.Number(1.6), incx.rhs()))
    elif cell_type == 1:
        incx.set_rhs(myokit.Multiply(myokit.Number(2.0), incx.rhs()))
    else:
        print 'Enter valid cell type for grandi model (Endo: 0, Epi: 1)'

    # Sarco/ endoplasmic reticulum Ca2+ pump current- GPB 50% in HF
    serca = m.get('caflux.J_serca')
    if cell_type == None:
        serca.set_rhs(myokit.Multiply(myokit.Number(0.5), serca.rhs()))
    elif cell_type == 0:
        serca.set_rhs(myokit.Multiply(myokit.Number(0.45), serca.rhs()))
    elif cell_type == 1:
        serca.set_rhs(myokit.Multiply(myokit.Number(0.75), serca.rhs()))
    else:
        print 'Enter valid cell type for grandi model (Endo: 0, Epi: 1)'


    # SR Ca2+ leak GPB 300% in HF
    leak = m.get('caflux.J_SRleak')
    leak.set_rhs(myokit.Multiply(myokit.Number(3.0), leak.rhs()))

    # EC50. -11% in HF
    ec50 = m.get('caflux.ec50SR')
    ec50.set_rhs(myokit.Multiply(myokit.Number(0.89), ec50.rhs()))

    return(m)

## Elshrif (2015) HF model, heterogenous across cell types
def Ord_HF_Elshrif(cell_type):
    m = myokit.load_model('ohara-2011.mmt')

    # Fast sodium channel
    ina = m.get('ina.INa')
    ina.set_rhs(myokit.Multiply(myokit.Number(0.37), ina.rhs()))

    # Late sodium current
    inal = m.get('inal.INaL')
    inal.set_rhs(myokit.Multiply(myokit.Number(1.93), inal.rhs()))

    # No time constant of inactivation change for INal in Elshrif model

    # Transient outward K+ current
    ito = m.get('ito.Ito')
    if cell_type == 1:
        ito.set_rhs(myokit.Multiply(myokit.Number(0.6), ito.rhs()))
    elif cell_type == 2:
        ito.set_rhs(myokit.Multiply(myokit.Number(0.62), ito.rhs()))
    elif cell_type == 0:
        ito.set_rhs(myokit.Multiply(myokit.Number(0.49), ito.rhs()))
    else:
        print 'Enter valid cell_type for Elshrif Ord HF model'

    # Inward rectifier K+ current
    ik1 = m.get('ik1.IK1')
    if cell_type == 0 or cell_type == 1:
        ik1.set_rhs(myokit.Multiply(myokit.Number(0.45), ik1.rhs()))
    elif cell_type == 2:
        ik1.set_rhs(myokit.Multiply(myokit.Number(0.47), ik1.rhs()))

    # Na+/ K+ pump current
    inak = m.get('inak.INaK')
    # All cell types -40%
    inak.set_rhs(myokit.Multiply(myokit.Number(0.6), inak.rhs()))

    # Sodium Calcium Exchanger
    incx = m.get('inaca.INaCa')
    # All cell types +131%
    incx.set_rhs(myokit.Multiply(myokit.Number(2.31), incx.rhs()))

    # SERCA pump
    serc = m.get('serca.Jup') # Could be . Jtr or others, not clear
    # Epi/ Endo
    if cell_type == 1 or cell_type == 0:
        serc.set_rhs(myokit.Multiply(myokit.Number(0.59), serc.rhs()))
    # Mid-myocardial
    elif cell_type == 2:
        serc.set_rhs(myokit.Multiply(myokit.Number(0.58), serc.rhs()))

    # IKr, rapid delayed rectifier potassium current
    ikr = m.get('ikr.IKr')
    # Epi
    if cell_type == 1:
        ikr.set_rhs(myokit.Multiply(myokit.Number(0.54), ikr.rhs()))
    # Endo
    elif cell_type == 0:
        ikr.set_rhs(myokit.Multiply(myokit.Number(0.73), ikr.rhs()))
    # Mid-myocardial cells- no change -> Don't need to multiply by value

    # IKs, slow delayed rectifier potassium current
    iks = m.get('iks.IKs')
    # Epi
    if cell_type == 1:
        iks.set_rhs(myokit.Multiply(myokit.Number(0.41), iks.rhs()))
    # Endo
    elif cell_type == 0:
        iks.set_rhs(myokit.Multiply(myokit.Number(0.42), iks.rhs()))
    # Mid-myocardial cells
    elif cell_type == 2:
        iks.set_rhs(myokit.Multiply(myokit.Number(0.5), iks.rhs()))

    return(m)


def GPB_HF_Moreno():
    m = myokit.load_model('grandi-2010.mmt')

    # Original grandi model does not contain late channel. Changes from Tenor et al.
    # Late Na+ current. +900%
    inal = m.get('inal.GNal') # GNal multiplied by both junc and sl
    inal.set_rhs(myokit.Multiply(myokit.Number(10.0), inal.rhs()))

    # Transient outward K+ current -36%
    ito = m.get('ito.ito')
    ito.set_rhs(myokit.Multiply(myokit.Number(0.64), ito.rhs()))

    # Inward rectifier K+ current. -25% in HF
    ik1 = m.get('ik1.I_k1')
    ik1.set_rhs(myokit.Multiply(myokit.Number(0.75), ik1.rhs()))

    # Na+/ K+ pump current- Between -42% and -10% in literature. Taken mean:-26%
    inak = m.get('inak.I_nak')
    inak.set_rhs(myokit.Multiply(myokit.Number(0.74), inak.rhs()))

    # Background Na+ current +1600% in moreno (2013) from rabbit HF data
    inab = m.get('inab.GNaB') # Then slow and junc = 0, as both multiplied by GNaB
    inab.set_rhs(myokit.Multiply(myokit.Number(17.0), inab.rhs()))

    # Sarco/ endoplasmic reticulum Ca2+ pump current: -36% in HF human patients
    serca = m.get('caflux.J_serca')
    serca.set_rhs(myokit.Multiply(myokit.Number(0.64), serca.rhs()))

    # SR Ca2+ leak. +350% from rabbit HF data
    leak = m.get('caflux.J_SRleak')
    leak.set_rhs(myokit.Multiply(myokit.Number(4.5), leak.rhs()))

    return(m)

def TT_HF_Lu():
    m = myokit.load_model('tentusscher-2006.mmt')

    # INa: -40%
    ina = m.get('ina.INa')
    ina.set_rhs(myokit.Multiply(myokit.Number(0.6), ina.rhs()))

    # Ito: -36%
    ito = m.get('ito.ITo')
    ito.set_rhs(myokit.Multiply(myokit.Number(0.64), ito.rhs()))

    # IK1: -20%
    ik1 = m.get('ik1.IK1')
    ik1.set_rhs(myokit.Multiply(myokit.Number(0.8), ik1.rhs()))

    # INaK: -42%
    inak = m.get('inak.INaK')
    inak.set_rhs(myokit.Multiply(myokit.Number(0.58), inak.rhs()))

    # ICab: +53%
    icab = m.get('icab.ICaB')
    icab.set_rhs(myokit.Multiply(myokit.Number(1.53), icab.rhs()))

    # INCX: +65%
    incx = m.get('inaca.INaCa')
    incx.set_rhs(myokit.Multiply(myokit.Number(1.65), incx.rhs()))


    # SR calcium pump current : -45%
    iup = m.get('calcium.i_up')
    iup.set_rhs(myokit.Multiply(myokit.Number(0.55), iup.rhs()))
    ## Can I call this serca?

    # ILeak: 0%
    ileak = m.get('calcium.i_leak')
    ileak.set_rhs(myokit.Multiply(myokit.Number(1.0), ileak.rhs()))


    #Jrel ??: -23%
    jrel = m.get('calcium.i_rel') # Not sure this is the right parameter
    jrel.set_rhs(myokit.Multiply(myokit.Number(0.77), jrel.rhs()))


    #Iks: -50%
    iks = m.get('iks.IKs')
    iks.set_rhs(myokit.Multiply(myokit.Number(0.5), iks.rhs()))

    return(m)

# Main function for testing
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
