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

# livetimes
livetimes = array([1497.44, 793.71, 562.04, 2790.1])

# energy-independent efficiency sys
# TODO: Update efficiencies for SK-IV
efftot = {"ando": [0.7975,0.56703,0.77969], "malaney": [0.7903, 0.53273, 0.766262]}
sys_eff = array([0.0254, 0.0404, 0.0253, 0.03])
sys_eff_sk4_ntag = 0.12
regionid = {"low": 0, "medium": 1, "high": 2}
pdfid = {"nue": 0, "numu": 1, "nc": 2, "mupi": 3, "rel": 4}
modelid = {"ando": 0}

# 3rd reduction efficiencies
effs_sk4 = loadtxt("efficiencies/efficiencies_sk4.txt")
effsk4 = interp1d(effs_sk4[:,0], effs_sk4[:,1], bounds_error = False, fill_value = (effs_sk4[0,1], effs_sk4[-1,1]))
effs_3rdred = [None, None, None, effsk4]

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
spaeff_sk2 = array([[16,0.762], [18,0.882], [24, 1.0]])
spaeff_sk3 = array([[16,0.818], [18,0.908], [24, 1.0]])

spaeff_sk4_ntag = array([[16,0.88],[18,0.88],[24,1.0]])

spaeff = [spaeff_sk1, spaeff_sk2, spaeff_sk3, spaeff_sk4]

# Solar efficiencies
soleff_sk4 = array([[16,0.736], [17,0.813], [18,0.865], [19, 0.965], [20, 1]])
soleff_sk1 = array([[16,0.738], [17,0.821], [18,0.878], [19, 0.965], [20, 1]])
soleff_sk2 = array([[17.5,0.738], [18.02,0.821], [19.08,0.878], [20.14, 0.965], [26,1]])
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

# Load background pdfs
# WARNING!!! 16 MeV threshold is hardcoded there!
bg_sk4_dir = "./pdf_bg_sk4"
cut_bins_ntag, cut_effs_ntag = get_spasolbins(4,ntag = True)
#cut_bins_ntag, cut_effs_ntag = [16, 18, 24, 90], [0.88, 0.88, 1.0]
bgs_sk4_ntag = [bg_sk4(i, cut_bins_ntag, cut_effs_ntag, bg_sk4_dir, 16., ntag = True) for i in range(4)]
cut_bins, cut_effs = get_spasolbins(4,ntag = False)
#cut_bins, cut_effs = [16, 17, 18, 19, 20, 24, 90], [0.373 * 0.736, 0.373 * 0.813, 0.705 * 0.865, 0.753 * 0.965, 0.866, 1.0]
bgs_sk4 = [bg_sk4(i, cut_bins, cut_effs, bg_sk4_dir, 16., ntag = False) for i in range(4)]

def get_eff(model, sknum, signal = None):
    if sknum < 4:
        eff = efftot[model][sknum - 1]
    else:
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

def load_signal_pdf_sk4(model, elow, ntag = False):
    print(f"Load model {model}")
    if ":" not in model:
        # These are discrete models
        flux = loadtxt(f"models/flux_cross_{model}.dat")
        en = arange(elow,90.1,0.1)
        spec = array([sns.ibd_spectrum_flux(ee, flux) for ee in en])
        return relic_sk4(en, spec, lambda z: seff_sk(z, elow, 4, ntag = ntag), elow)
    else:
        raise ValueError("Parameterized models not implemented yet!")

def pdf(energy, sknum, model, elow, pdfid, region, signal = None, ntag = False):
    if sknum < 4: return likes.pdf(energy, sknum, modelid[model], elow, pdfid, region)
    elif sknum == 4:
        if pdfid == 4: return signal.pdf(energy,region) # to do (for specific srn models)
        elif pdfid in range(4) and ntag: return bgs_sk4_ntag[pdfid].pdf(energy, region)
        elif pdfid in range(4) and not ntag: return bgs_sk4[pdfid].pdf(energy, region)
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
def systematics_atm(energies_med, energies_high, sknum, model, elow):
    sigmas = arange(-1, 3.5, 0.5)
    # Normalization and correction factors for nue CC
    norm0 = quad(lambda en: pdf(en, sknum, model, elow, pdfid["nue"], regionid["medium"]), elow, 90)[0]
    norm1 = quad(lambda en: en * pdf(en, sknum, model, elow, pdfid["nue"], regionid["medium"]), elow, 90)[0]
    normnue = 1./(1 + 0.5 * sigmas * (norm1/norm0 - 16)/74)
    nuefact = 1 + 0.5 * sigmas[newaxis,:] * (energies_med[:,newaxis] - 16)/74
    nuefact *= normnue
    # Correction factors for NC
    normncmed = quad(lambda en: pdf(en, sknum, model, elow, pdfid["nc"], regionid["medium"]), elow, 90)[0]
    normnchigh = quad(lambda en: pdf(en, sknum, model, elow, pdfid["nc"], regionid["high"]), elow, 90)[0]
    print(normncmed, normnchigh)
    ncfact_med = 1 + sigmas
    ncfact_high = 1 - sigmas * normncmed/normnchigh
    ncfact_high = where(ncfact_high < 0, 0, ncfact_high)
    print(ncfact_high, normnue)
    # make systematics tensors
    sysmatrix_med = ones((len(energies_med), 5, len(sigmas), len(sigmas)))
    #sysmatrix_med[:,pdfid["nue"],:,:] = nuefact[...,newaxis] 
    #sysmatrix_med[:,pdfid["nc"],:,:] = broadcast_to(ncfact_med, (len(energies_med), len(sigmas), len(sigmas)))

    sysmatrix_high = ones((len(energies_high), 5, len(sigmas), len(sigmas)))
    #sysmatrix_high[:,pdfid["nc"],:,:] = broadcast_to(ncfact_high, (len(energies_high), len(sigmas), len(sigmas)))
    return sysmatrix_med, sysmatrix_high

# Asymmetric gaussian for atm systematics weighting
def asym_gaussian():
    return array([0.1643, 0.2517, 0.2636, 0.1888, 0.09240, 0.03092, 0.007076, 0.001107, 0.0001184])

# Get pdfs for different energies, regions, types
# Output is an array of pdf values for each energy and each type of signal/bg 
def get_pdfmatrix(energies, sknum, region, model, elow, signal = None, ntag = False):
    if sknum in [1, 2, 3]:
        return array(likes.get_pdfs(energies, sknum, regionid[region], modelid[model], elow))
    else:
        p = [[pdf(e, sknum, model, elow, i, regionid[region], signal = signal, ntag = ntag) for i in range(5)] for e in energies]
        return array(p)

# Likelihood without systematics
def get_like_nosys(ncce, nnc, nmupi, nrelic, pdfs_low, pdfs_med, pdfs_high):
    ntot = len(pdfs_low) + len(pdfs_high) + len(pdfs_med)
    nccmu = ntot - ncce - nnc - nmupi - nrelic
    if ncce < 0 or nccmu < 0: return -1e10
    nevents = array([ncce[0], nccmu] + [nnc, nmupi] + [nrelic])
    totlike = log((nevents * pdfs_med).sum(axis = 1)).sum() - nevents.sum()
    print(nevents.sum(), totlike + nevents.sum(), totlike)
    #print ncce, nccmu, pdfs_med[30], nevents,log((nevents * pdfs_med).sum(axis = 1))[30]
    return totlike

# Likelihood with systematics
def get_like(nbackgrounds, nrelic, pdfs_distorted_low, pdfs_distorted_med, pdfs_distorted_high, sysmatrix_med, sysmatrix_high):
    if nbackgrounds.min() < 0: return -1e10
    wgauss = asym_gaussian()
    nevents = array(list(nbackgrounds) + [nrelic])
    totlike = (log(einsum("j,ijkl", nevents, pdfs_distorted_high)).sum(axis = 0) 
               + log(dot(nevents,pdfs_distorted_low.T)).sum(axis = 0) 
               + log(einsum("j,ijkl", nevents, pdfs_distorted_med)).sum(axis = 0) 
               - nrelic - nbackgrounds.sum())
    totmax = totlike.max()
    likenew = log((exp(totlike - totmax) * wgauss[:, newaxis] * wgauss[newaxis,:]).sum()) + totmax
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
def maxlike(sknum, model, elow, ntag = False, rmin = -5, rmax = 100, rstep = 0.1):
    likemax = -1e10
    # Get signal spectrum
    signal = None
    effsignal = 1.0
    if sknum== 4:
        signal = load_signal_pdf_sk4(model,elow,ntag)
        effsignal = get_eff(model,sknum,signal)
    else:
        effsignal = efftot[model][sknum - 1]
    print(f"Efficiency is {effsignal}")

    # Get pdfs
    samplow, sampmed, samphigh = (low[sknum - 1], med[sknum - 1], high[sknum - 1]) if not ntag else (low_ntag, med_ntag, high_ntag)
    samplow = random_sample(len(samplow)) * 74 + 16
    pdfs_high = get_pdfmatrix(samphigh, sknum, "high", model, elow, signal,ntag)
    pdfs_med = get_pdfmatrix(sampmed, sknum, "medium", model, elow, signal,ntag)
    pdfs_low = get_pdfmatrix(samplow, sknum, "low", model, elow, signal,ntag)
    #print(pdfs_high)

    # Get systematic error matrices
    sysmatrix_med, sysmatrix_high = systematics_atm(sampmed, samphigh, sknum, model, elow)

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
    #print("90% c.l. {} /cm^2/s {} MeV".format(limit * flux(model), elow))
    lpos = results_sys[:, -1].argmax()
    print("nu-e events {}".format(results_sys[lpos, 1]))
    print("nu-mu events {}".format(results_sys[lpos, 2]))
    print("NC elastic events {}".format(results_sys[lpos, 3]))
    print("mu/pi events {}".format(results_sys[lpos, 4]))
    print("Max likelihood {}".format(results_sys[lpos, -1]))
    print("Numbers of events: {} {} {}".format(len(samplow), len(sampmed), len(samphigh)))
    plotfit(results_sys[lpos,1], results_sys[lpos,2], results_sys[lpos,3], results_sys[lpos,4], results_sys[lpos,5], model, sknum, elow, signal = signal, ntag = ntag)

def plotfit(nnue, nnumu, nnc, nmupi, nrelic, model, sknum, elow, signal = None, ntag = False):
    #fig, ax = plt.subplots(1,3,sharey = True)
    samplow, sampmed, samphigh = (low[sknum - 1], med[sknum - 1], high[sknum - 1]) if not ntag else (low_ntag, med_ntag, high_ntag)
    en= arange(elow, 90, 0.1)
    nuecc = nnue * array([pdf(ee, sknum, model, elow, pdfid["nue"], 1) for ee in en])
    numucc = nnumu * array([pdf(ee, sknum, model, elow, pdfid["numu"], 1) for ee in en])
    nc = nnc * array([pdf(ee, sknum, model, elow, pdfid["nc"], 1) for ee in en])
    mupi = nmupi * array([pdf(ee, sknum, model, elow, pdfid["mupi"], 1) for ee in en])
    plt.plot(en, nuecc, label = "nue CC")
    plt.plot(en, nc, label = "NC")
    plt.plot(en, numucc, label = "numu CC")
    plt.plot(en, mupi, label = "mupi")
    plt.plot(en, mupi + nc + numucc + nuecc, label = "All")
    h = histogram(sampmed, bins = arange(16,90,7.4), weights = 1./7.4 * ones(len(sampmed)))
    plt.errorbar(arange(16+7.4/2,90,7.4)[:-1], h[0], xerr = 7.4/2, yerr = sqrt(h[0]), fmt = '.', color = 'black', label = "Data")
    plt.legend()
    plt.xlabel("Ep (MeV)")
    plt.ylabel("Number of events after cuts")
    # High
    plt.figure()
    en= arange(elow, 90, 0.1)
    nuecc = nnue * array([pdf(ee, sknum, model, elow, pdfid["nue"], 2) for ee in en])
    numucc = nnumu * array([pdf(ee, sknum, model, elow, pdfid["numu"], 2) for ee in en])
    nc = nnc * array([pdf(ee, sknum, model, elow, pdfid["nc"], 2) for ee in en])
    mupi = nmupi * array([pdf(ee, sknum, model, elow, pdfid["mupi"], 2) for ee in en])
    plt.plot(en, nuecc, label = "nue CC")
    plt.plot(en, nc, label = "NC")
    plt.plot(en, numucc, label = "numu CC")
    plt.plot(en, mupi, label = "mupi")
    plt.plot(en, mupi + nc + numucc + nuecc, label = "All")
    h = histogram(samphigh, bins = arange(16,90,7.4), weights = 1./7.4 * ones(len(samphigh)))
    plt.errorbar(arange(16+7.4/2,90,7.4)[:-1], h[0], xerr = 7.4/2, yerr = sqrt(h[0]), fmt = '.', color = 'black', label = "Data")
    plt.legend()
    plt.xlabel("Ep (MeV)")
    plt.ylabel("Number of events after cuts")
    # Low
    plt.figure()
    en= arange(elow, 90, 0.1)
    nuecc = nnue * array([pdf(ee, sknum, model, elow, pdfid["nue"], 0) for ee in en])
    numucc = nnumu * array([pdf(ee, sknum, model, elow, pdfid["numu"], 0) for ee in en])
    nc = nnc * array([pdf(ee, sknum, model, elow, pdfid["nc"], 0) for ee in en])
    mupi = nmupi * array([pdf(ee, sknum, model, elow, pdfid["mupi"], 0) for ee in en])
    plt.plot(en, nuecc, label = "nue CC")
    plt.plot(en, nc, label = "NC")
    plt.plot(en, numucc, label = "numu CC")
    plt.plot(en, mupi, label = "mupi")
    plt.plot(en, mupi + nc + numucc + nuecc, label = "All")
    h = histogram(samplow, bins = arange(16,90,7.4), weights = 1./7.4 * ones(len(samplow)))
    plt.errorbar(arange(16+7.4/2,90,7.4)[:-1], h[0], xerr = 7.4/2, yerr = sqrt(h[0]), fmt = '.', color = 'black', label = "Data")
    plt.legend()
    plt.xlabel("Ep (MeV)")
    plt.ylabel("Number of events after cuts")

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
