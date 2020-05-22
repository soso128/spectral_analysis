from __future__ import division
from numpy import *
from scipy.interpolate import interp1d
from scipy.integrate import quad
from scipy.optimize import fmin
from pdf_sk4 import bg_sk4, relic_sk4
from sys import path
path.append("spectrum_generator/")
import snspectrum as sns
import likes

# livetimes
livetimes = array([1497.44, 793.71, 562.04, 2790.1])

# energy-independent efficiency sys
# TODO: Update efficiencies for SK-IV
efftot = {"ando": [0.7975,0.56703,0.77969], "malaney": [0.7903, 0.53273, 0.766262]}
sys_eff = array([0.0254, 0.0404, 0.0253, 0.03])
regionid = {"low": 0, "medium": 1, "high": 2}
pdfid = {"nue": 0, "numu": 1, "nc": 2, "mupi": 3, "rel": 4}
modelid = {"ando": 0}

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

# Load background pdfs
bg_sk4_dir = "./pdf_bg_sk4"
cut_bins, cut_effs = [16, 18, 24, 90], [0.88, 0.88, 1.0]
bgs_sk4 = [bg_sk4(i, cut_bins, cut_effs, bg_sk4_dir, 16., ntag = True) for i in range(4)]

def get_eff(model, sknum, signal = None):
    if sknum < 4:
        eff = efftot[model][sknum - 1]
    else:
        eff = signal.overall_efficiency()
    return eff

# Signal efficiencies for SK-IV
def seff_sk4(en, elow):
    if en < elow: return 0
    if en < 18: return 0.231 * 0.88
    if en < 24: return 0.298 * 0.88
    if en < 30: return 0.283
    if en < 40: return 0.263
    if en < 60: return 0.269
    if en < 90: return 0.283
    return 0

def load_signal_pdf_sk4(model, elow):
    print(f"Load model {model}")
    if ":" not in model:
        # These are discrete models
        flux = loadtxt(f"models/flux_cross_{model}.dat")
        en = arange(elow,90.1,0.1)
        spec = array([sns.ibd_spectrum_flux(ee, flux) for ee in en])
        return relic_sk4(en, spec, lambda z: seff_sk4(z, elow), elow)
    else:
        raise ValueError("Parameterized models not implemented yet!")

def pdf(energy, sknum, model, elow, pdfid, region, signal = None):
    if sknum < 4: return likes.pdf(energy, sknum, model, elow, pdfid, region)
    elif sknum == 4:
        if pdfid == 4: return signal.pdf(energy,region) # to do (for specific srn models)
        elif pdfid in range(4): return bgs_sk4[pdfid].pdf(energy, region)
        else: raise ValueError("Invalid pdfid")
    else: raise ValueError("Invalid sknum")

# Samples for SK I-IV
def load_sample(sknum):
    low = loadtxt("sk{}/samplelow.txt".format(int(sknum)))[:, 1]
    med = loadtxt("sk{}/samplemed.txt".format(int(sknum)))[:, 1]
    high = loadtxt("sk{}/samplehigh.txt".format(int(sknum)))[:, 1]
    return low,med,high

low1, med1, high1 = load_sample(1)
low2, med2, high2 = load_sample(2)
low3, med3, high3 = load_sample(3)
low4, med4, high4 = load_sample(4)
#low4, med4, high4 = load_sample(4)
low = [low1, low2, low3, low4]
med = [med1, med2, med3, med4]
high = [high1, high2, high3, high4]

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
    ncfact_med = 1 + sigmas
    ncfact_high = 1 - sigmas * normncmed/normnchigh
    print(ncfact_high, normnue)
    # make systematics tensors
    sysmatrix_med = ones((len(energies_med), 5, len(sigmas), len(sigmas)))
    sysmatrix_med[:,pdfid["nue"],:,:] = nuefact[...,newaxis] 
    sysmatrix_med[:,pdfid["nc"],:,:] = broadcast_to(ncfact_med, (len(energies_med), len(sigmas), len(sigmas)))

    sysmatrix_high = ones((len(energies_high), 5, len(sigmas), len(sigmas)))
    sysmatrix_high[:,pdfid["nc"],:,:] = broadcast_to(ncfact_high, (len(energies_high), len(sigmas), len(sigmas)))
    return sysmatrix_med, sysmatrix_high

# Asymmetric gaussian for atm systematics weighting
def asym_gaussian():
    return array([0.1643, 0.2517, 0.2636, 0.1888, 0.09240, 0.03092, 0.007076, 0.001107, 0.0001184])

# Get pdfs for different energies, regions, types
# Output is an array of pdf values for each energy and each type of signal/bg 
def get_pdfmatrix(energies, sknum, region, model, elow, signal = None):
    if sknum in [1, 2, 3]:
        return array(likes.get_pdfs(energies, sknum, regionid[region], modelid[model], elow))
    else:
        p = [[pdf(e, sknum, model, elow, i, regionid[region], signal) for i in range(5)] for e in energies]
        return array(p)

# Likelihood without systematics
def get_like_nosys(ncce, nnc, nmupi, nrelic, sknum, model, elow, pdfs_low, pdfs_med, pdfs_high):
    ntot = len(pdfs_low) + len(pdfs_high) + len(pdfs_med)
    nccmu = ntot - ncce - nnc - nmupi - nrelic
    if ncce < 0 or nccmu < 0: return -1e10
    nevents = array([ncce[0], nccmu] + [nnc, nmupi] + [nrelic])
    totlike = log((nevents * pdfs_med).sum(axis = 1)).sum() - nevents.sum()
    #print ncce, nccmu, pdfs_med[30], nevents,log((nevents * pdfs_med).sum(axis = 1))[30]
    return totlike

# Likelihood with systematics
def get_like(nbackgrounds, nrelic, sknum, model, elow, pdfs_distorted_low, pdfs_distorted_med, pdfs_distorted_high, sysmatrix_med, sysmatrix_high):
    if nbackgrounds.min() < 0: return -1e10
    wgauss = asym_gaussian()
    nevents = array(list(nbackgrounds) + [nrelic])
    #print(where(einsum("j,ijkl", nevents, pdfs_distorted_high) < 0))
    totlike = (log(einsum("j,ijkl", nevents, pdfs_distorted_high)).sum(axis = 0) 
               + log(dot(nevents,pdfs_distorted_low.T)).sum(axis = 0) 
               + log(einsum("j,ijkl", nevents, pdfs_distorted_med)).sum(axis = 0) 
               - nrelic - nbackgrounds.sum())
    totmax = totlike.max()
    likenew = log((exp(totlike - totmax) * wgauss[:, newaxis] * wgauss[newaxis,:]).sum()) + totmax
    return likenew

def getmaxlike(nrelic, nback_ini, sknum, model, elow, pdfs_low, pdfs_med, pdfs_high, sys = 0, sysmatrix_med = None, sysmatrix_high  = None):
    funclike = None
    if sys:
        funclike = lambda nback: -get_like(nback, nrelic, sknum, model, elow, pdfs_low, pdfs_med, pdfs_high, sysmatrix_med, sysmatrix_high)
        maxlike = fmin(funclike, nback_ini, full_output = True, disp = 0)
        return append(maxlike[0], array([nrelic, -maxlike[1]]))
    else:
        funclike = lambda nback: -get_like_nosys(nback, nback_ini[1], nback_ini[2], nrelic, sknum, model, elow, pdfs_low, pdfs_med, pdfs_high)
        maxlike = fmin(funclike, nback_ini[0], full_output = True, disp = 0)
        ntot = len(pdfs_low) + len(pdfs_high) + len(pdfs_med)
        nccmu = ntot - maxlike[0] - nback_ini[1] - nback_ini[2] - nrelic
        return concatenate([array([maxlike[0][0], nccmu]), nback_ini[1:], array([nrelic, -maxlike[1]])])

# Main maximum likelihood function
# sknum = 1,2,3 (SK phase)
# model = SRN model (right now, only "ando")
# elow = energy threshold (here, 16MeV)
# rmin, rmax, rstep = range and step of numbers of relic events for likelihood maximization
def maxlike(sknum, model, elow, rmin = -5, rmax = 100, rstep = 0.1):
    likemax = -1e10
    # Get signal spectrum
    signal = None
    effsignal = 1.0
    if sknum== 4:
        signal = load_signal_pdf_sk4(model,elow)
        effsignal = get_eff(model,sknum,signal)
    else:
        effsignal = efftot[model][sknum - 1]

    # Get pdfs
    pdfs_high = get_pdfmatrix(high[sknum - 1], sknum, "high", model, elow, signal)
    pdfs_med = get_pdfmatrix(med[sknum - 1], sknum, "medium", model, elow, signal)
    pdfs_low = get_pdfmatrix(low[sknum - 1], sknum, "low", model, elow, signal)

    # Get systematic error matrices
    sysmatrix_med, sysmatrix_high = systematics_atm(med[sknum-1], high[sknum - 1], sknum, model, elow)

    # Distort pdfs
    pdfs_syshigh = pdfs_high[...,newaxis,newaxis] * sysmatrix_high
    pdfs_sysmed = pdfs_med[...,newaxis,newaxis] * sysmatrix_med

    # Set backgrounds (preliminary estimates)
    nue, numu, nc, mupi, relic, loglik = initialize(sknum, model, elow, pdfs_low, pdfs_med, pdfs_high, rmin)

    # Main maximization loop
    likedata = []
    rmin = 0
    for i,rel in enumerate(arange(rmin, rmax, rstep)):
        likedata.append(getmaxlike(rel, array([nue, numu, nc, mupi]), sknum, model, elow, pdfs_low, pdfs_sysmed, pdfs_syshigh, sys = 1, sysmatrix_med = sysmatrix_med, sysmatrix_high = sysmatrix_high))
        # Update initial values
        [nue, numu, nc, mupi] = list(likedata[-1][:4])
        if i % 100 == 0: print("Step {}/1000, like = {}".format(i, likedata[-1][-1]))
        
    results = column_stack((arange(rmin, rmax, rstep), likedata))
    results = results[results[:, 0] >= 0]

    # Systematic efficiency error correction + limits
    lmax2, best2, errplus2, errminus2, limit2 = analyse(results)
    results_sys = applysys(results, sknum, model,effsignal)
    lmax, best, errplus, errminus, limit = analyse(results_sys, final = True)

    # Save and display results
    savetxt("fit_sk{}.txt".format(sknum), column_stack((results, results_sys[:, -1])))
    print("SK-{}".format(sknum), "Best fit:")
    print("{} +{} -{} relic evts/yr".format(best, errplus, errminus))
    print("{} +{} -{} relic evts".format(best2, errplus2, errminus2))
    print("90% c.l. relic event rate: {} ev/yr {}", limit, limit2)
    #print("90% c.l. {} /cm^2/s {} MeV".format(limit * flux(model), elow))
    lpos = results_sys[:, -1].argmax()
    print("nu-e events {}".format(results_sys[lpos, 1]))
    print("nu-mu events {}".format(results_sys[lpos, 2]))
    print("NC elastic events {}".format(results_sys[lpos, 3]))
    print("mu/pi events {}".format(results_sys[lpos, 4]))
    print("Max likelihood {}".format(results_sys[lpos, -1]))
    print("Numbers of events: {} {} {}".format(len(low[sknum - 1]), len(med[sknum - 1]), len(high[sknum - 1])))

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

def applysys(likes,sknum,model,eff):
    # Get gaussian pdf for systematics
    lower = max(eff - 6*sys_eff[sknum - 1], 0)
    upper = min(eff + 6*sys_eff[sknum - 1], 1)
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
def initialize(sknum, model, ethr, low, med, high, rel):
    nevlow = len(low)
    nevhigh = len(high)
    nevmed = len(med)
    nback = nevlow + nevhigh + nevmed - rel

    # Estimate mupi and nc backgrounds using sidebands
    # Fractions are taken from MC
    mupi = nevlow * mupi_rescale_low[sknum - 1]
    nc = (nevhigh - mupi_rescale_high[sknum - 1] * mupi) * nc_rescale[sknum - 1]

    # Maximize likelihoods over backgrounds in signal region
    likemax = getmaxlike(rel, array([nback/5.,nc,mupi]), sknum, model, ethr, low, med, high)

    return likemax

#maxlike(3, "ando", 16)
