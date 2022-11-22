from numpy import *
from scipy.integrate import quad
from scipy.interpolate import interp1d

# Energy scale and resolution systematics for the SRN analysis
res_sys = 0.03
scale_sys = 0.02

# Energy resolution functions in different phases of SK
def energy_sigma(ee, sknum):
    if sknum == 1:
        return 0.2468 + 0.1492 * sqrt(ee) + 0.0690 * ee
    if sknum == 2:
        return -0.123 + 0.376 * sqrt(ee) + 0.0349 * ee
    if sknum == 3:
        return 0.0536 + 0.5200 * sqrt(ee) + 0.0458 * ee
    if sknum >= 4:
        return -0.05525 + 0.3162 * sqrt(ee) + 0.04572 * ee
        #return -0.290 + 0.434 * sqrt(ee) + 0.0320 * ee

def get_res_smear(elow, eup, fspec, sknum):
    x = arange(elow, eup, 0.1)
    specnew = array([quad(lambda y: fspec(y) * exp(-(xx - y)**2/(4 * res_sys * energy_sigma(y, sknum)**2)), max(xx - 10 *sqrt(res_sys) * energy_sigma(xx, sknum), 0), min(eup + 10, xx + 10 *sqrt(res_sys) * energy_sigma(xx, sknum)))[0] for xx in x])
    normnew = array([quad(lambda y: exp(-(xx - y)**2/(4 * res_sys * energy_sigma(y, sknum)**2)), max(xx - 10 *sqrt(res_sys) * energy_sigma(xx, sknum), 0), min(eup + 10, xx + 10 *sqrt(res_sys) * energy_sigma(xx, sknum)))[0] for xx in x])
    specnew /= normnew
    specold = array([fspec(xx) for xx in x])
    xpol = specold > 0.001 * specold.max()
    spec_ratio = (specnew/specold)[xpol]
    pol = polyfit(x[xpol], spec_ratio, 3)
    return pol

def get_scale_smear(elow, eup, fspec, sknum):
    x = arange(elow, eup, 0.1)
    specnew = array([fspec(xx * (1 + scale_sys)) for xx in x])
    specold = array([fspec(xx) for xx in x])
    normnew = specnew.sum()/specold.sum()
    specnew /= normnew
    xpol = specold > 0.001 * specold.max()
    spec_ratio = (specnew/specold)[xpol]
    pol = polyfit(x[xpol], spec_ratio, 3)
    return pol

# Distorsion functions for SK phases >= 4
def get_distorsion_functions(elow_0n, ehigh, elow_1n, backgrounds, signal, sknum):
    scale_pdict = {}
    res_pdict = {}
    # Signal (consider only medium angle region)
    elow = [elow_0n, elow_1n]
    for ntag in [0,1]:
        fspec = lambda x: signal.pdf(x, 1, ntag)
        pscale = get_scale_smear(elow[ntag], ehigh, fspec, sknum)
        pres = get_res_smear(elow[ntag], ehigh, fspec, sknum)
        scale_pdict[f"relic{ntag}"] = pscale
        res_pdict[f"relic{ntag}"] = pres
    # Backgrounds
    bg_names = ["ccnue", "decaye", "nc", "mupi"]
    region_names = ["low", "med", "high"] 
    for bg,bgname in zip(backgrounds, bg_names):
        for region,regname in zip([0,1,2], region_names):
            for ntag in [0,1]:
                fspec = lambda x: bg.pdf(x, region, ntag)
                pscale = get_scale_smear(elow[ntag], ehigh, fspec, sknum)
                pres = get_res_smear(elow[ntag], ehigh, fspec, sknum)
                scale_pdict[f"{bgname}{regname}{ntag}"] = pscale
                res_pdict[f"{bgname}{regname}{ntag}"] = pres
    return scale_pdict, res_pdict
