from __future__ import division
from numpy import *
from numpy.random import random_sample
from scipy.interpolate import interp1d
from scipy.integrate import quad
from scipy.optimize import fmin
from pdf_sk4 import bg_sk4, relic_sk4
from sys import path
path.append("spectrum_generator/")
import snspectrum as sns
import likes
import matplotlib.pyplot as plt
from ast import literal_eval
import seaborn
import matplotlib as mpl

mpl.rcParams['font.size'] = 15

# livetimes
livetimes = array([1497.44, 793.71, 562.04, 2790.1])

# energy-independent efficiency sys
# TODO: Update efficiencies for SK-IV
efftot = {"lma": [0.7975,0.56703,0.77969], "malaney": [0.7903, 0.53273, 0.766262], "faild": [0.7997,0.5845,0.7831] , "woosley": [0.7876,0.5524,0.7764], "ksw" :[0.7908, 0.54527, 0.77371]}
fluxfac = {"lma": 0.535115, "malan": 0.5845, "ksw": 0.488413, "faild": 0.431, "woosley": 0.47447}
sys_eff = array([0.0254, 0.0404, 0.0253, 0.03])
sys_eff_sk4_ntag = 0.12
regionid = {"low": 0, "medium": 1, "high": 2}
pdfid = {"nue": 0, "numu": 1, "nc": 2, "mupi": 3, "rel": 4}
modelid = {"lma": 0, "faild": -3, "malaney": -1, "ksw": -2, "woosley": -4}

# 3rd reduction efficiencies
effs_sk4 = loadtxt("efficiencies/efficiencies_sk4.txt")
effsk4 = interp1d(effs_sk4[:,0], effs_sk4[:,1], bounds_error = False, fill_value = (effs_sk4[0,1], effs_sk4[-1,1]))
effs_sk3 = loadtxt("efficiencies/efficiencies_sk3.txt")
effsk3 = interp1d(effs_sk3[:,0], effs_sk3[:,1], bounds_error = False, fill_value = (effs_sk3[0,1], effs_sk3[-1,1]))
effs_sk2 = loadtxt("efficiencies/efficiencies_sk2.txt")
effsk2 = interp1d(effs_sk2[:,0], effs_sk2[:,1], bounds_error = False, fill_value = (effs_sk2[0,1], effs_sk2[-1,1]))
effs_sk1 = loadtxt("efficiencies/efficiencies_sk1.txt")
effsk1 = interp1d(effs_sk1[:,0], effs_sk1[:,1], bounds_error = False, fill_value = (effs_sk1[0,1], effs_sk1[-1,1]))
effs_3rdred = [effsk1, effsk2, effsk3, effsk4]

# Spallation + 3rd red efficiency curves as a function of energy for SK I-IV
#eff_sk1 = loadtxt("input/eff_sk1.dat")
#eff_sk2 = loadtxt("input/eff_sk2.dat")
#eff_sk3 = loadtxt("input/eff_sk3.dat")
#eff_sk4 = loadtxt("input/eff_sk4.dat")
#effsk1 = interp1d(eff_sk1[:, 0], eff_sk1[:, 1])
#effsk2 = interp1d(eff_sk2[:, 0], eff_sk2[:, 1])
#effsk3 = interp1d(eff_sk3[:, 0], eff_sk3[:, 1])
#effsk4 = interp1d(eff_sk4[:, 0], eff_sk4[:, 1])

#effs = [effsk1, effsk2, effsk3, effsk4]

# Scalings between Cherenkov angle regions (from MC)
mupi_rescale_low = [1.367, 1.75, 1.34, 1.34] # mupi from low to medium
mupi_rescale_high = [0.12777, 0.1, 0.13, 0.13] # mupi from low to high
nc_rescale = [1.16313, 1.42, 1.14, 1.14] # NC from high to medium

# Spallation efficiencies
spaeff_sk4 = array([[16, 0.373], [17, 0.373], [18, 0.705], [19, 0.753], [20, 0.866], [24, 1.0]])
spaeff_sk1 = array([[16,0.818], [18,0.908], [24, 1.0]])
spaeff_sk2 = array([[17.5,0.762], [20,0.882], [26, 1.0]])
spaeff_sk3 = array([[16,0.818], [18,0.908], [24, 1.0]])

spaeff_sk4_ntag = array([[16,0.88],[18,0.88],[24,1.0]])

spaeff = [spaeff_sk1, spaeff_sk2, spaeff_sk3, spaeff_sk4]

# Solar efficiencies
soleff_sk4 = array([[16,0.736], [17,0.813], [18,0.865], [19, 0.965], [20, 1]])
soleff_sk1 = array([[16,0.738], [17,0.821], [18,0.878], [19, 0.965], [20, 1]])
soleff_sk2 = array([[17.5,0.738], [18.02,0.821], [19.08,0.878], [20.14, 0.965], [21.2,1]])
soleff_sk3 = array([[16,0.738], [17,0.821], [18,0.878], [19, 0.965], [20, 1]])

soleff = [soleff_sk1, soleff_sk2, soleff_sk3, soleff_sk4]

# Get efficiency steps for background pdfs
def get_spasolbins(sknum, ntag = False):
    bins = []
    effs = []
    if sknum < 4 or not ntag:
        spa = spaeff[sknum- 1]
        sol = soleff[sknum - 1]
        spabins = spa[:, 0]
        solbins = sol[:, 0]
        effspa = lambda x: 1.0 if x < spa[0,0] else spa[(spa[:, 0] <= x), -1][-1]
        effsol = lambda x: 1.0 if x < sol[0,0] else sol[(sol[:, 0] <= x), -1][-1]
        bins = array(list(set(sorted(append(spabins, solbins)))))
        effs = vectorize(effspa)(bins) * vectorize(effsol)(bins)
        bins = list(bins) + [90.]
    elif sknum == 4 and ntag:
        bins = list(spaeff_sk4_ntag[:, 0]) + [90.]
        effs = list(spaeff_sk4_ntag[:, 1])
    else:
        raise ValueError(f"No search for SK-{sknum} and ntag={ntag}")
    return bins, effs


def get_eff(model, sknum, signal = None):
    #if sknum < 4:
        #eff = efftot[model][sknum - 1]
    #else:
        #eff = signal.overall_efficiency()
    eff = signal.overall_efficiency()
    return eff

# Signal efficiencies (solar + spallation + 3rd red + ntag)
# No solar cut with ntag
def seff_sk(en, elow, sknum, ntag = False):
    if en < elow: return 0
    eff3 = effs_3rdred[sknum - 1](en)
    if ntag:
        if sknum < 4:
            raise ValueError(f"SK-{sknum} does not have ntag")
        spa = spaeff_sk4_ntag
        effspa = lambda x: 1.0 if x < spa[0,0] else spa[(spa[:, 0] <= x), -1][-1]
        eff0 = eff3 * effspa(en)
        if en < 18: return eff0 * .231
        if en < 24: return eff0 * .298
        if en < 30: return eff0 * .283
        if en < 40: return eff0 * .263
        if en < 60: return eff0 * .269
        if en < 90: return eff0 * .283
    else:
        spa = spaeff[sknum - 1]
        sol = soleff[sknum - 1]
        effspa = lambda x: 1.0 if x < spa[0,0] else spa[(spa[:, 0] <= x), -1][-1]
        effsol = lambda x: 1.0 if x < sol[0,0] else sol[(sol[:, 0] <= x), -1][-1]
        return eff3 * effspa(en) * effsol(en)
    return 0

def load_signal_pdf(sknum,model, elow, ntag = False):
    print(f"Load model {model}")
    if ":" not in model:
        # These are discrete models
        flux = loadtxt(f"models/flux_cross_{model}.dat")
        en = arange(elow,90.1,0.1)
        spec = array([sns.ibd_spectrum_flux(ee, flux) for ee in en])
        ee,spec = sns.smear_ibd_spectrum(en,column_stack((en,spec)),sknum)
        return relic_sk4(en, spec, lambda z: seff_sk(z, elow, sknum, ntag = ntag), elow),fluxfac[model]
    else:
        print(f"Loading {model}")
        imf = sns.imfs['salpeter']
        csfr = sns.csfr_fiducial
        mod = literal_eval(model)
        en = arange(elow,90.1,0.1)
        spec = array([sns.ibd_spectrum(ee, mod, imf, csfr) for ee in en])
        print(spec)
        ee,spec = sns.smear_ibd_spectrum(en,column_stack((en,spec)),sknum)
        totrate = spec[ee>elow].sum() * (ee[1] - ee[0])
        return relic_sk4(en, spec, lambda z: seff_sk(z, elow, sknum, ntag = ntag), elow),totrate

def pdf(energy, sknum, model, elow, pdfid, region, backgrounds = None, signal = None):
    if sknum < 4: 
        if pdfid == 4: return likes.pdf(energy, sknum, modelid[model], elow, 4, region)
        #if pdfid == 4: return signal.pdf(energy,region)
        if pdfid in range(4): return likes.pdf(energy, sknum, 0, elow, pdfid, region)
    elif sknum == 4:
        if pdfid == 4: return signal.pdf(energy,region) # to do (for specific srn models)
        elif pdfid in range(4): return backgrounds[pdfid].pdf(energy, region)
        else: raise ValueError("Invalid pdfid")
    else: raise ValueError("Invalid sknum")

# Samples for SK I-IV
def load_sample(sknum, ntag = False):
    low = loadtxt("sk{}/samplelow.txt".format(int(sknum)))[:, 1]
    med = loadtxt("sk{}/samplemed.txt".format(int(sknum)))[:, 1]
    high = loadtxt("sk{}/samplehigh.txt".format(int(sknum)))[:, 1]
    if ntag:
        low = loadtxt("sk{}/ntag/samplelow.txt".format(int(sknum)))[:, 1]
        med = loadtxt("sk{}/ntag/samplemed.txt".format(int(sknum)))[:, 1]
        high = loadtxt("sk{}/ntag/samplehigh.txt".format(int(sknum)))[:, 1]
    return low,med,high

low1, med1, high1 = load_sample(1)
low2, med2, high2 = load_sample(2)
low3, med3, high3 = load_sample(3)
low4, med4, high4 = load_sample(4)
#low4, med4, high4 = load_sample(4)
low = [low1, low2, low3, low4]
med = [med1, med2, med3, med4]
high = [high1, high2, high3, high4]

# ntag samples for SK-IV only
low_ntag, med_ntag, high_ntag = load_sample(4, ntag = True)

# Compute distorsion functions due to systematics (for atmospheric spectra)
def systematics_atm(energies_med, energies_high, sknum, model, elow, backgrounds = None):
    sigmas = arange(-1, 3.5, 0.5)
    sigmas2 = arange(-1, 3.5, 0.5)
    # Normalization and correction factors for nue CC
    norm0 = quad(lambda en: pdf(en, sknum, model, elow, pdfid["nue"], regionid["medium"], backgrounds), elow, 90)[0]
    norm1 = quad(lambda en: en * pdf(en, sknum, model, elow, pdfid["nue"], regionid["medium"], backgrounds), elow, 90)[0]
    normnue = 1./(1 + 0.5 * sigmas * (norm1/norm0 - 16)/74)
    nuefact = 1 + 0.5 * sigmas[newaxis,:] * (energies_med[:,newaxis] - 16)/74
    nuefact *= normnue
    # Correction factors for NC
    normncmed = quad(lambda en: pdf(en, sknum, model, elow, pdfid["nc"], regionid["medium"], backgrounds), elow, 90)[0]
    normnchigh = quad(lambda en: pdf(en, sknum, model, elow, pdfid["nc"], regionid["high"], backgrounds), elow, 90)[0]
    if sknum == 4: sigmas2 = arange(-2,2,0.5)/2.
    ncfact_med = 1 + sigmas2
    ncfact_high = 1 - sigmas2 * normncmed/normnchigh
    #ncfact_high = where(ncfact_high < 0, 0, ncfact_high)
    print(ncfact_high, normnue)
    # make systematics tensors
    sysmatrix_med = ones((len(energies_med), 5, len(sigmas), len(sigmas2)))
    sysmatrix_med[:,pdfid["nue"],:,:] = nuefact[...,newaxis] 
    sysmatrix_med[:,pdfid["nc"],:,:] = broadcast_to(ncfact_med, (len(energies_med), len(sigmas), len(sigmas2)))

    sysmatrix_high = ones((len(energies_high), 5, len(sigmas), len(sigmas2)))
    sysmatrix_high[:,pdfid["nc"],:,:] = broadcast_to(ncfact_high, (len(energies_high), len(sigmas), len(sigmas2)))
    return sysmatrix_med, sysmatrix_high

# Asymmetric gaussian for atm systematics weighting
def asym_gaussian():
    return array([0.1643, 0.2517, 0.2636, 0.1888, 0.09240, 0.03092, 0.007076, 0.001107, 0.0001184])

# Get pdfs for different energies, regions, types
# Output is an array of pdf values for each energy and each type of signal/bg 
def get_pdfmatrix(energies, sknum, region, model, elow, signal = None, backgrounds = None):
    #if sknum in [1, 2, 3]:
        #return array(likes.get_pdfs(energies, sknum, regionid[region], modelid[model], elow))
    #else:
    p = [[pdf(e, sknum, model, elow, i, regionid[region], signal = signal, backgrounds = backgrounds) for i in range(5)] for e in energies]
    return array(p)

# Likelihood without systematics
# Only used for initialization of background numbers
def get_like_nosys(ncce, nnc, nmupi, nrelic, pdfs_low, pdfs_med, pdfs_high):
    ntot = len(pdfs_low) + len(pdfs_high) + len(pdfs_med)
    nccmu = ntot - ncce - nnc - nmupi - nrelic
    if ncce < 0 or nccmu < 0: return -1e10
    nevents = array([ncce[0], nccmu] + [nnc, nmupi] + [nrelic])
    totlike = log((nevents * pdfs_med).sum(axis = 1)).sum() - nevents.sum()
    #print ncce, nccmu, pdfs_med[30], nevents,log((nevents * pdfs_med).sum(axis = 1))[30]
    return totlike

# Likelihood with systematics
def get_like(nbackgrounds, nrelic, pdfs_distorted_low, pdfs_distorted_med, pdfs_distorted_high, sysmatrix_med, sysmatrix_high):
    if nbackgrounds.min() < 0: return -1e10
    wgauss = asym_gaussian() 
    sknum = 3 if sysmatrix_med.shape[-2] == sysmatrix_med.shape[-1] else 4
    wgauss2 = asym_gaussian() if sknum < 4 else 0.20997 * exp(-arange(-2,2,0.5)**2/2.)
    nevents = array(list(nbackgrounds) + [nrelic])
    totlike = (log(einsum("j,ijkl", nevents, pdfs_distorted_high)).sum(axis = 0) 
               + log(dot(nevents,pdfs_distorted_low.T)).sum(axis = 0) 
               + log(einsum("j,ijkl", nevents, pdfs_distorted_med)).sum(axis = 0) 
               - nrelic - nbackgrounds.sum()) # maybe double counting?
    totmax = totlike.max()
    likenew = log((exp(totlike - totmax) * wgauss[:, newaxis] * wgauss2[newaxis,:]).sum()) + totmax
    return likenew

def getmaxlike(nrelic, nback_ini, pdfs_low, pdfs_med, pdfs_high, sys = 0, sysmatrix_med = None, sysmatrix_high  = None):
    funclike = None
    if sys:
        funclike = lambda nback: -get_like(nback, nrelic, pdfs_low, pdfs_med, pdfs_high, sysmatrix_med, sysmatrix_high)
        maxlike = fmin(funclike, nback_ini, full_output = True, disp = 0)
        return append(maxlike[0], array([nrelic, -maxlike[1]]))
    else:
        funclike = lambda nback: -get_like_nosys(nback, nback_ini[1], nback_ini[2], nrelic, pdfs_low, pdfs_med, pdfs_high)
        maxlike = fmin(funclike, nback_ini[0], full_output = True, disp = 0)
        ntot = len(pdfs_low) + len(pdfs_high) + len(pdfs_med)
        nccmu = ntot - maxlike[0] - nback_ini[1] - nback_ini[2] - nrelic
        return concatenate([array([maxlike[0][0], nccmu]), nback_ini[1:], array([nrelic, -maxlike[1]])])

# Main maximum likelihood function
# sknum = 1,2,3 (SK phase)
# model = SRN model (right now, only "ando")
# elow = energy threshold (here, 16MeV)
# rmin, rmax, rstep = range and step of numbers of relic events for likelihood maximization
def maxlike(sknum, model, elow, ntag = False, rmin = -5, rmax = 100, rstep = 0.1, quiet = True):
    likemax = -1e10
    # Get signal spectrum
    signal = None
    effsignal = 1.0
    totrate = 0
    if sknum== 4:
        signal,totrate = load_signal_pdf(sknum,model,elow,ntag)
        effsignal = get_eff(model,sknum,signal)
    else:
        signal,totrate = load_signal_pdf(sknum,model,elow,ntag)
        effsignal = get_eff(model,sknum,signal)
        effsignal = efftot[model][sknum - 1]
    print(f"Efficiency is {effsignal}")

    # Get background spectra
    bgs_sk4 = None
    bg_sk4_dir = "./pdf_bg_sk4"
    if sknum == 4 and ntag:
        # Load background pdfs
        # WARNING!!! 16 MeV threshold is hardcoded there!
        cut_bins_ntag, cut_effs_ntag = get_spasolbins(4,ntag = True)
        #cut_bins_ntag, cut_effs_ntag = [16, 18, 24, 90], [0.88, 0.88, 1.0]
        bgs_sk4 = [bg_sk4(i, cut_bins_ntag, cut_effs_ntag, bg_sk4_dir, elow, ntag = True) for i in range(4)]
    else:
        cut_bins, cut_effs = get_spasolbins(4,ntag = False)
        #cut_bins, cut_effs = [16, 17, 18, 19, 20, 24, 90], [0.373 * 0.736, 0.373 * 0.813, 0.705 * 0.865, 0.753 * 0.965, 0.866, 1.0]
        bgs_sk4 = [bg_sk4(i, cut_bins, cut_effs, bg_sk4_dir, elow, ntag = False) for i in range(4)]
        

    # Get pdfs
    samplow, sampmed, samphigh = (low[sknum - 1], med[sknum - 1], high[sknum - 1]) if not ntag else (low_ntag, med_ntag, high_ntag)
    pdfs_high = get_pdfmatrix(samphigh, sknum, "high", model, elow, signal = signal, backgrounds = bgs_sk4)
    pdfs_med = get_pdfmatrix(sampmed, sknum, "medium", model, elow, signal= signal, backgrounds = bgs_sk4)
    pdfs_low = get_pdfmatrix(samplow, sknum, "low", model, elow, signal = signal, backgrounds = bgs_sk4)
    #print(pdfs_high)

    # Get systematic error matrices
    sysmatrix_med, sysmatrix_high = systematics_atm(sampmed, samphigh, sknum, model, elow, backgrounds = bgs_sk4)

    # Distort pdfs
    pdfs_syshigh = pdfs_high[...,newaxis,newaxis] * sysmatrix_high
    pdfs_sysmed = pdfs_med[...,newaxis,newaxis] * sysmatrix_med

    # Set backgrounds (preliminary estimates)
    nue, numu, nc, mupi, relic, loglik = initialize(sknum, pdfs_low, pdfs_med, pdfs_high, rmin)

    # Main maximization loop
    likedata = []
    rmin = 0
    for i,rel in enumerate(arange(rmin, rmax, rstep)):
        likedata.append(getmaxlike(rel, array([nue, numu, nc, mupi]), pdfs_low, pdfs_sysmed, pdfs_syshigh, sys = 1, sysmatrix_med = sysmatrix_med, sysmatrix_high = sysmatrix_high))
        # Update initial values
        [nue, numu, nc, mupi] = list(likedata[-1][:4])
        if i % 100 == 0: print("Step {}/1000, like = {}".format(i, likedata[-1][-1]))
        
    results = column_stack((arange(rmin, rmax, rstep), likedata))
    results = results[results[:, 0] >= 0]

    # Systematic efficiency error correction + limits
    lmax2, best2, errplus2, errminus2, limit2 = analyse(results)
    results_sys = applysys(results, sknum, effsignal,ntag)
    lmax, best, errplus, errminus, limit = analyse(results_sys, final = True)

    # Save and display results
    savetxt("fit_sk{}.txt".format(sknum), column_stack((results, results_sys[:, -1])))
    print("SK-{}".format(sknum), "Best fit:")
    print("{} +{} -{} relic evts/yr".format(best, errplus, errminus))
    print("{} +{} -{} relic evts".format(best2, errplus2, errminus2))
    print("90% c.l. relic event rate: {} ev/yr {}".format(limit, limit2))
    if ':' not in model:
        print("90% c.l. {} /cm^2/s {} MeV".format(limit * fluxfac[model], elow))
    lpos = results_sys[:, -1].argmax()
    print("nu-e events {}".format(results_sys[lpos, 1]))
    print("nu-mu events {}".format(results_sys[lpos, 2]))
    print("NC elastic events {}".format(results_sys[lpos, 3]))
    print("mu/pi events {}".format(results_sys[lpos, 4]))
    print("Max likelihood {}".format(results_sys[lpos, -1]))
    print("Numbers of events: {} {} {}".format(len(samplow), len(sampmed), len(samphigh)))
    if not quiet:
        plotfit(results_sys[lpos,1], results_sys[lpos,2], results_sys[lpos,3], results_sys[lpos,4], results_sys[lpos,5], model, sknum, elow, signal = signal, background = bgs_sk4, ntag = ntag)
    if ':' not in model: return limit*fluxfac[model],fluxfac[model],results_sys
    return limit,totrate,results_sys

def fullike(model, elow, rmin = -5, rmax = 100, rstep = 0.1, quiet = False):
    like1 = maxlike(1, model, elow, rmin = rmin, rmax = rmax, rstep = rstep)
    like2 = maxlike(2, model, elow, rmin = rmin, rmax = rmax, rstep = rstep)
    like3 = maxlike(3, model, elow, rmin = rmin, rmax = rmax, rstep = rstep)
    like4 = maxlike(4, model, elow, rmin = rmin, rmax = rmax, rstep = rstep)
    like4ntag = maxlike(4, model, elow, ntag = True, rmin = rmin, rmax = rmax, rstep = rstep)
    res = combine([like1[-1],like2[-1],like3[-1], like4[-1], like4ntag[-1]])
    if not quiet:
        plt.figure()
        plt.xlabel("SRN events/year")
        plt.ylabel("Likelihood")
        x = like1[-1][:, 0]/2
        plt.plot(x, like1[-1][:, -1] - like1[-1][:,-1].max(), label = "SK-I", alpha = 0.5)
        plt.plot(x, like2[-1][:, -1] - like2[-1][:,-1].max(),'--', label = "SK-II", alpha = 0.5)
        plt.plot(x, like3[-1][:, -1] - like3[-1][:,-1].max(),'-.', label = "SK-III", alpha = 0.5)
        plt.plot(x, like4ntag[-1][:, -1] - like4ntag[-1][:,-1].max(), '-.', label = "SK-IV (ntag)", linewidth=2)
        plt.plot(x, like4[-1][:, -1] - like4[-1][:,-1].max(), '--', label = "SK-IV",linewidth=2)
        likesum = like1[-1][:, -1] + like2[-1][:, -1] + like3[-1][:, -1] + like4[-1][:, -1] + like4ntag[-1][:, -1]
        likesum -= likesum.max()
        plt.plot(x, likesum, label = "Combined", color = 'black')
        plt.plot([0,20], -0.5 * ones(2), 'r')
        plt.xlim(0,20)
        plt.ylim(-2,0.2)
        plt.grid()
        plt.legend()
    if ":" not in model: res = list(res) + [like1[1] * res[-1]]
    print(f"Limits are {like1[0]} {like2[0]} {like3[0]} {like4[0]} {like4ntag[0]}")
    return list(res) + [like1[1]]

def fullike_sk4(model, elow, rmin = -5, rmax = 100, rstep = 0.1, quiet = False):
    like4 = maxlike(4, model, elow, rmin = rmin, rmax = rmax, rstep = rstep)
    like4ntag = maxlike(4, model, elow, ntag = True, rmin = rmin, rmax = rmax, rstep = rstep)
    res = combine([like4[-1], like4ntag[-1]])
    if not quiet:
        plt.figure()
        plt.xlabel("SRN events/year")
        plt.ylabel("Likelihood")
        x = like4[-1][:, 0]/2
        plt.plot(x, like4ntag[-1][:, -1] - like4ntag[-1][:,-1].max(),'--', label = "SK-IV (ntag)", linewidth=2)
        plt.plot(x, like4[-1][:, -1] - like4[-1][:,-1].max(),':', label = "SK-IV",linewidth=2)
        likesum = like4[-1][:, -1] + like4ntag[-1][:, -1]
        likesum -= likesum.max()
        plt.plot(x, likesum, label = "Combined", color = 'black')
        plt.plot([0,20], -0.5 * ones(2), 'r')
        plt.xlim(0,20)
        plt.ylim(-2,0.2)
        plt.grid()
        plt.legend()
    if ":" not in model: res = list(res) + [like4[1] * res[-1]]
    print(f"Limits are {like4[0]} {like4ntag[0]}")
    return list(res) + [like4[1]]

def combine(results):
    liketot = results[0][:, -1] - results[0][:, -1].max()
    for r in results[1:]:
        liketot += r[:, -1] - r[:, -1].max()
    return analyse(column_stack((results[0][:,:-1], liketot)), final = True)

def plotregion(region, data, nnue, nnumu, nnc, nmupi, nrelic, model, sknum, elow, ax, signal = None, background = None):
    #plt.figure()
    step = 2
    en= arange(elow, 90, 0.1)
    nuecc = nnue * array([pdf(ee, sknum, model, elow, pdfid["nue"], region, background) for ee in en])
    numucc = nnumu * array([pdf(ee, sknum, model, elow, pdfid["numu"], region, background) for ee in en])
    nc = nnc * array([pdf(ee, sknum, model, elow, pdfid["nc"], region, background) for ee in en])
    mupi = nmupi * array([pdf(ee, sknum, model, elow, pdfid["mupi"], region, background) for ee in en])
    relic = nrelic * array([pdf(ee, sknum, model, elow, 4, region, signal = signal) for ee in en])
    ax.plot(en, step*nuecc, label = r"$\nu_e$ CC")
    ax.plot(en, step*nc, label = "NC")
    ax.plot(en, step*numucc, label = r"$\nu_\mu$ CC")
    ax.plot(en, step*mupi, label = r"$\mu/\pi$")
    ax.plot(en, step*relic, label = "relic")
    ax.plot(en, step*(mupi + nc + numucc + nuecc), label = "all background")
    h = histogram(data, bins = arange(elow,90,step))
    x = h[1][1:] - 0.5 * (h[1][1:] - h[1][:-1])
    ax.errorbar(x, h[0], xerr = step/2, yerr = sqrt(h[0]), fmt = '.', color = 'black')
    if ax.is_first_col():
        ax.legend(loc = 'upper left', prop={'size': 12})
        ax.set(ylabel = "Number of events after cuts")
    ax.set(xlabel = "Ep (MeV)")

def plotfit(nnue, nnumu, nnc, nmupi, nrelic, model, sknum, elow, signal = None, background = None, ntag = False):
    fig, (ax1, ax2, ax3) = plt.subplots(1,3,sharey = True)
    samplow, sampmed, samphigh = (low[sknum - 1], med[sknum - 1], high[sknum - 1]) if not ntag else (low_ntag, med_ntag, high_ntag)
    plotregion(1, sampmed, nnue, nnumu, nnc, nmupi, nrelic, model, sknum, elow, ax2, signal, background)
    plotregion(2, samphigh, nnue, nnumu, nnc, nmupi, nrelic, model, sknum, elow, ax3, signal, background)
    plotregion(0, samplow, nnue, nnumu, nnc, nmupi, nrelic, model, sknum, elow, ax1, signal, background)
    plt.subplots_adjust(wspace = 0)

def analyse(likes, final = False):
    lmax = likes[:, -1].max()
    bestpos = likes[:, -1].argmax()
    fact = 0.5 if final else 1
    rel = likes[:, -2] * fact
    best = rel[bestpos]
    norm = exp(likes[:, -1] - lmax).sum()
    errminus = best - rel[searchsorted(likes[:bestpos, -1], lmax - 0.5)]
    errplus = rel[len(likes) - searchsorted(likes[bestpos:, -1][::-1], lmax - 0.5)] - best
    l90 = rel[searchsorted(exp(likes[:, -1] - lmax).cumsum(), 0.9 * norm)] 
    return lmax, best, errplus, errminus, l90

def applysys(likes,sknum,eff,ntag = False):
    # Get gaussian pdf for systematics
    print(f"Signal efficiency is {eff}")
    if sknum < 4 and ntag:
        raise ValueError(f"No ntag with SK-{sknum}")
    syseff = sys_eff[sknum - 1] if not ntag else sys_eff_sk4_ntag
    lower = max(eff * (1 - 6*syseff), 1e-10)
    upper = min(eff * (1 + 6*syseff), 0.999)
    step = (upper - lower)/1000.
    epsrange = arange(lower, upper+step, step)
    if len(epsrange) > 1001: epsrange = epsrange[:-1]
    pgaus = exp(-0.5 * (epsrange - eff)**2/(sys_eff[sknum - 1]*eff)**2)
    pgaus /= pgaus.sum()
    # Convolution (Simpson integration)
    lmax = likes[:, -1].max()
    flikes = interp1d(likes[:, -2], exp(likes[:, -1] - lmax), bounds_error = False, fill_value = 0)
    #print(flikes(10))
    #print(exp(likes[:, -1] - lmax))
    rates = arange(0,100,0.1)
    lconv = flikes(rates[:, newaxis] * epsrange[newaxis, :] * livetimes[sknum - 1]/365.25* 0.5)
    simpsoncoeff = array([step/3.] + list((1 + (arange(1,1000)%2))*2./3 * step) + [step/3.])
    ltot = (lconv * (pgaus * epsrange * simpsoncoeff)).sum(axis = 1)
    likenew = log(ltot) + lmax
    return column_stack((likes[:, :-1], likenew))

# Initialization
def initialize(sknum, low, med, high, rel):
    nevlow = len(low)
    nevhigh = len(high)
    nevmed = len(med)
    nback = nevlow + nevhigh + nevmed - rel

    # Estimate mupi and nc backgrounds using sidebands
    # Fractions are taken from MC
    mupi = nevlow * mupi_rescale_low[sknum - 1]
    nc = (nevhigh - mupi_rescale_high[sknum - 1] * mupi) * nc_rescale[sknum - 1]

    # Maximize likelihoods over backgrounds in signal region
    likemax = getmaxlike(rel, array([nback/5.,nc,mupi]), low, med, high)
    print(mupi, nc)
    print(likemax)
    print(nevlow, nevhigh, nevmed)

    return likemax

#maxlike(3, "ando", 16)
