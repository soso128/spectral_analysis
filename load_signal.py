from numpy import *
from scipy.interpolate import interp1d
from scipy.integrate import quad
from sys import path


from fit_parameters import *
path.append("spectrum_generator/")
import snspectrum as sns
import get_flux_from_mc as gtmc
from pdf_sk4 import  relic_sk

def load_signal_pdf(sknum, model, elow, ehigh, elow_1n, ineff_scale=1.0):
    ''' Load relic signal model pdfs, either from a discrete (named) model
    or a parametric model. '''

    def seff_sk(en, ntag=False):
        ''' Signal efficiencies (solar + spallation + 3rd red + ntag)
        (No solar cut with ntag) '''
        if (not ntag and en < elow) or (ntag and en < elow_1n):
            return 0
        eff3 = effs_3rdred[sknum - 1](en)
        if ntag:
            if sknum < 4:
                raise ValueError(f"SK-{sknum} does not have ntag")
            spa = spaeff_sk4
            effspa = lambda x: 1.0 if x < spa[0,0] else spa[(spa[:, 0] <= x), -1][-1]
            tot_seff = eff3 * effspa(en)
            return tot_seff + (1-tot_seff)*(1-ineff_scale)
        else:
            spa = spaeff[sknum - 1]
            sol = soleff[sknum - 1]
            effspa = lambda x: 1.0 if x < spa[0,0] else spa[(spa[:, 0] <= x), -1][-1]
            effsol = lambda x: 1.0 if x < sol[0,0] else sol[(sol[:, 0] <= x), -1][-1]
            tot_seff = eff3 * effspa(en) * effsol(en)
            return tot_seff + (1-tot_seff)*(1-ineff_scale)
        return 0

    def relic(en, spec):
        eff_func_nontag = lambda z: seff_sk(z, ntag=False)
        if sknum < 4:
            return relic_sk(sknum, en, spec, s_ch_frac, eff_func_nontag,
                            elow=elow, ehigh=ehigh)
        else:
            # print("ntag_effs:", ntag_effs * ntag_eff_ps)
            # print("ntag_bgs:", ntag_bgs * ntag_bg_ps)
            # print("ntag_eff_ps:", ntag_eff_ps)
            # print("ntag_bg_ps:", ntag_bg_ps)
            eff_func_ntag = lambda z: seff_sk(z, ntag=True)
            return relic_sk(sknum, en, spec, s_ch_frac, eff_func_nontag,
                            eff_func_ntag, elow=elow, ehigh=ehigh,
                            elow_n=elow_1n, ntag_ebins=ntag_ebins,
                            ntag_effs=ntag_effs, ntag_bgs=ntag_bgs,
                            ntag_eff_ps=ntag_eff_ps, ntag_bg_ps=ntag_bg_ps)

    def spec_flux(flux, low, high, smear=True):
        """ Return spectrum given neutrino flux,
        and positron energy bounds """
        en = arange(low, high + 0.1, 0.1)
        spec0 = array([sns.ibd_spectrum_flux(ee, flux) for ee in en])
        spec, en = spec0[~isnan(spec0)], en[~isnan(spec0)] # remove undefined
        if smear:
            specfunc = gtmc.smeared_spec_from_mc(en, spec, sknum, high)
            spec = specfunc(en)
            #_, spec = sns.smear_ibd_spectrum(en, column_stack((en, spec)), sknum)
        return en, spec

    def spec_params(mod, imf, csfr, low, high, smear=True):
        """ Return spectrum given parametrization,
        and positron energy bounds """
        en = arange(low, high + 0.1, 0.1)
        spec0 = array([sns.ibd_spectrum(ee, mod, imf, csfr) for ee in en])
        spec, en = spec0[~isnan(spec0)], en[~isnan(spec0)] # remove undefined
        if smear:
            specfunc = gtmc.smeared_spec_from_mc(en, spec, sknum, high)
            spec = specfunc(en)
            #_, spec = sns.smear_ibd_spectrum(en, column_stack((en, spec)), sknum)
        return en, spec

    def tot_rate(en_full, spec_full, low, high, ntag=None, nfracs=None):
        en = en_full[(en_full >= low) & (en_full <= high)]
        spec = spec_full[(en_full >= low) & (en_full <= high)]
        if ntag is not None:
            totrate = 0
            for e, sp in zip(en, spec):
                ebin = digitize(e, ntag_ebins) - 1
                nfrac = nfracs[ntag][ebin]
                totrate += nfrac * sp
        else:
            totrate = spec.sum()
        return totrate * (en[1] - en[0])


    print(f"Load model {model}")
    if ":" not in model:
        # These are discrete models
        flux = loadtxt(f"models/flux_cross_{model}.dat")  # Enu, flux
        fflux = interp1d(flux[:, 0], flux[:, 1], bounds_error=False, fill_value=0)
        totflux, _ = quad(fflux, 17.3, 100)
        en, spec_full = spec_flux(flux, 0, 100)
    else:
        print(f"Loading {model}")
        imf = sns.imfs['salpeter']
        csfr = sns.csfr_fiducial
        mod = literal_eval(model)
        print(mod, imf, csfr)
        totflux, _ = quad(sns.snflux, 17.3, 100, args=(mod, imf, csfr))
        en, spec_full = spec_params(mod, imf, csfr, 0, 100)
    rel = relic(en, spec_full)
    totrate = tot_rate(en, spec_full, 16, 90)
    flux_fac = totflux / totrate
    return rel, flux_fac, totrate, totflux