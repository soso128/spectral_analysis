from numpy import *
import snspectrum as sn
from scipy.interpolate import interp1d
from scipy.integrate import quad
#import ROOT as r

# Load MC ROOT files for SN models with blackbody neutrino spectra
# 2 arguments here: SK phase and neutrino temperature
def loadmc_root(tnu, sknum):
    data = None
    if sknum < 4:
        ch = r.TChain("h1")
        #tname = f"{int(tnu)}"
        tname="{}".format(int(tnu))
        if int(tnu) < tnu - 0.1: tname += ".5"
        ch.Add("../mc_root/relic{}mev_sk{}.root".format(tname, sknum))
        data = array([[ev.e,ev.mce,ev.x,ev.y,ev.z,ev.good,1] for ev in ch])
    else:
        ch = r.TChain("data")
        ch.Add("../mc_root/mc_sk4_formatted.root")
        data = array([[ev.bsenergy, ev.positron_energy, ev.x, ev.y, ev.z, ev.good, ev.dirks] for ev in ch])
        data[:, -1] = data[:, -2]**2 - data[:, -1]**2
    print("Data loaded")
    return data


def savemc(tnu, sknum):
    data = loadmc_root(tnu, sknum)
    if sknum < 4:
        tname="{}".format(int(tnu))
        if int(tnu) < tnu - 0.1: tname += ".5"
        save("../mc_npy/relic{}mev_sk{}.npy".format(tname, sknum), data)
    else:
        save("../mc_npy/mc_sk4_formatted.npy", data)

def loadmc(tnu, sknum):
    data = None
    if sknum < 4:
        tname="{}".format(int(tnu))
        if tnu > float(tname) + 0.1: tname += ".5"
        data = load("mc_npy/relic{}mev_sk{}.npy".format(tname, sknum), data)
    else:
        data = load("mc_npy/mc_sk4_formatted.npy", data)
    return data

#for sknum in range(1,4):
    #for tnu in arange(2,8,0.5):
        #savemc(tnu, sknum)

#savemc(1,4)


# First (noise) reduction
# Assume the columns have been reordered for each model file
def firstred(f, sknum):
    r = sqrt(f[:, 2]**2 + f[:, 3]**2)
    fv = (r < 1490) & (abs(f[:, 4]) < 1610)
    # OvaQ is not computed in SK-I
    ovaq = (f[:, 0] > -1) if sknum == 1 else (f[:, 6] > 0.25)
    return f[fv & ovaq & (f[:, 5] > 0.5) & (f[:, 0] < 9000)]

def smeared_spec_from_mc(energ, rates, sknum, eup = 90):
    #energ = arange(1,eup,0.1)
    ## Get rate for model
    #fflux = loadtxt(f"../models/flux_cross_{model}.dat")
    #fluxfunc = interp1d(fflux[:, 0], fflux[:, 1], bounds_error = False, fill_value = 0)
    #rates = array([sn.ibd_spectrum_flux(en,fflux) for en in energ])
    rate_func = interp1d(energ, rates, bounds_error = False, fill_value = 0)
    f = None
    blackbody_func = None
    if sknum < 4:
        # Sum all blackbody MC simulations
        blackbody_rates = array([array([sn.ibd_spectrum(en, {"type": 0, "tnu": tnu, "lumi": 3e53 * 1e3/0.0016}, sn.imfs['salpeter'], sn.csfr_fiducial) for en in energ]) for tnu in arange(3,7,0.5)])
        norm_rates = rates/rates.sum()
        blackbody_rates = blackbody_rates.T/blackbody_rates.sum(axis = 1)
        logdiff = ((log(norm_rates)[:,newaxis] - log(blackbody_rates))**2)[norm_rates > 0].sum(axis = 0)
        model_blackbody = blackbody_rates[:, logdiff.argmin()]
        #Load MC for blackbody model and reweight
        tnu = logdiff.argmin()*0.5 + 3
        #print(tnu)
        if (sknum == 1 and int(tnu + 0.1) >= 7): tnu = 6.5 # Smoother limits
        #if (int(tnu + 0.1) < 4.5): tnu = 3.5 # Smoother limits
        f = loadmc(tnu, sknum)
        blackbody_func = interp1d(energ, model_blackbody, bounds_error = False, fill_value = 0)
        #files = [loadmc(tnu, sknum) for tnu in arange(3,7,0.5)]
        #for ff in files: print(ff.shape)
        #f = concatenate(files, axis = 0)
        #blackbody_func = interp1d(energ, (blackbody_rates * array([len(ff) for ff in files])).sum(axis = 1)/len(f), bounds_error = False, fill_value = 0)
        #print(f.shape)
    else:
        f = loadmc(1, 4)
        blackbody_func = lambda x: ones(len(x))
    f1 = firstred(f,sknum)
    # Get rate for positron energies > 15 MeV by reweighting MC
    smeared_rate_high, _ = histogram(f1[:, 0], weights = rate_func(f1[:, 1])/blackbody_func(f1[:, 1]), bins = arange(0, eup, 1.0))
    true_rate_high, _ = histogram(f1[:, 1], weights = rate_func(f1[:, 1])/blackbody_func(f1[:, 1]), bins = arange(15, eup, 1.0))
    smeared_rate_high *= rates[energ > 15].sum()/true_rate_high.sum() * (energ[1] - energ[0])
    # Get rate at lower energies from resolution function
    smeared_rate_gauss = sn.smear_ibd_spectrum(energ, column_stack((energ,rates)),sknum)[1]
    # Make interpolation function
    rate_high_func = interp1d(arange(0.5,eup-1.0,1.0), smeared_rate_high, bounds_error = False, fill_value = 0)
    rate_low_func = interp1d(energ, smeared_rate_gauss, bounds_error = False, fill_value = 0)
    smeared_rate_func = lambda x: where(x < 15, rate_low_func(x), rate_high_func(x))
    return smeared_rate_func
