'''
Wrapper classes for relic and SK-IV backgrounds.
'''
from pickle import load as loadpick
from numpy import digitize, array, ones, exp
from scipy import interpolate
from scipy.integrate import quad


class bg_sk4:
    ''' Usage:
    nc = bg_sk4(2, cut_edges, efficiencies)
    nc.pdf(energy, region)
    '''
    def __init__(self, ev_type, cut_bins, cut_effs, cut_bins_n, cut_effs_n,
                 pdf_dir, elow, elow_n=None, ehigh=90.0, ntag_scale=None):
        ''' 0 <= ev_type <= 3 for nue, numu, nc, or mupi background.
        cut_bins should be a list of N energy bin edges.
        cut_effs should be a list of N-1 cut efficiencies.
        '''
        if elow < 10.0:
            raise ValueError("Can't go lower than 10 MeV!")
        if elow_n is None:
            elow_n = elow
        if elow_n < 10.0:
            raise ValueError("Can't go lower than 10 MeV!")
        self.elows = [elow, elow_n]
        self.ehigh = ehigh

        bins, effs = self._adjust_effs(cut_bins, cut_effs, elow, self.ehigh)
        bins_n, effs_n = self._adjust_effs(cut_bins_n, cut_effs_n, elow_n, self.ehigh)
        self.cut_bins = [bins, bins_n]
        self.cut_effs = [effs, effs_n]
        self.ev_type = ev_type
        self.pdf_dir = pdf_dir

        if ntag_scale is None:
            self.ntag_scale = ones(2)
        else:
            self.ntag_scale = ntag_scale[ev_type]
        self.params = self._load_param_array(self.ev_type)
        self.norm0 = self._get_norm0() # Normalization of source pdf to 1
        self.norm = self._get_norm() # Normalization after cuts to 1

    def _adjust_effs(self, bins, effs, elow, ehigh):
        newbins, neweffs = bins, effs
        if elow < bins[0]:
            newbins, neweffs = [elow] + list(bins), [1.0] + list(effs)
        if bins[-1] < ehigh:
            newbins, neweffs = list(bins) + [ehigh], list(effs) + [1.0]
        return newbins, neweffs

    def _load_param_array(self, ev_type):

        def _read_tck(file):
            with open(file, 'rb') as fl:
                tck = loadpick(fl, encoding = "bytes")
            return tck

        ev_labels = ["cc_e", "cc_mu", "nc", "mupi"]
        label = ev_labels[ev_type]
        files = ['%s/%s%d.p' % (self.pdf_dir, label, i) for i in range(3)]
        files_n = ['%s/%s%d_ntag.p' % (self.pdf_dir, label, i) for i in range(3)]

        params = [_read_tck(fl) for fl in files]
        params_n = [_read_tck(fl) for fl in files_n]
        return [params, params_n]

    def _pdf_before_cuts_raw(self, energy, region, ntag):
        ''' Unnormed pdf '''
        p = interpolate.splev(energy, self.params[ntag][region], der=0)
        return p.clip(min=1e-10) # Enforce positivity

    def _get_norm0(self):
        ''' Get pdf normalization before cuts '''
        area = 0.0
        for ntag in [0, 1]:
            for region in range(3):
                area += quad(self._pdf_before_cuts_raw, self.elows[ntag],
                             self.ehigh, args=(region, ntag))[0]
        return 1.0 / area

    def _check_valid_energy(self, energy, ntag):
        if energy < self.elows[ntag] or energy > self.ehigh:
            raise ValueError("Energy (%0.2f) outside range (%0.2f-%0.2f)"
                             % (energy, self.elows[ntag], self.ehigh))

    def pdf_before_cuts(self, energy, region, ntag):
        ''' Properly normalized pdf, before cuts '''
        self._check_valid_energy(energy, ntag)
        return self._pdf_before_cuts_raw(energy, region, ntag) * self.norm0

    def _get_norm(self):
        ''' Get pdf normalization after cuts '''
        area = 0.0
        for ntag in [0, 1]:
            for region in range(3):
                for i, cut_eff in enumerate(self.cut_effs[ntag]):
                    e_lo, e_hi = self.cut_bins[ntag][i], self.cut_bins[ntag][i+1]
                    if e_hi < self.elows[ntag]:
                        continue
                    elif e_lo < self.elows[ntag]:
                        e_lo = self.elows[ntag]
                    if e_lo > self.ehigh:
                        continue
                    elif e_hi > self.ehigh:
                        e_hi = self.ehigh
                    a0, _ = quad(self.pdf_before_cuts, e_lo,
                                 e_hi, args=(region, ntag))
                    area += a0 * cut_eff * self.ntag_scale[ntag]
        return 1.0 / area

    def pdf(self, energy, region, ntag):
        ''' Properly normalized pdf, after cuts '''
        try:
            self._check_valid_energy(energy, ntag)
        except ValueError:
            return 1e-10

        ebin = digitize(energy, self.cut_bins[ntag])
        if ebin < 1 or ebin >= len(self.cut_bins[ntag]):
            return 1e-10

        cut_eff = self.cut_effs[ntag][ebin-1]
        p0 = self.pdf_before_cuts(energy, region, ntag)
        return self.norm * cut_eff * p0 * self.ntag_scale[int(ntag)]

# Models spallation backgrounds
class spall_sk():
    ''' Usage (for 16 MeV threshold):
    spall = spall_sk4(cut_edges, efficiencies, efficiencies_3rd, 16) 
    spall.pdf(energy, region)
    '''
    def __init__(self, cut_bins, cut_effs, efficiency_func, sknum, elow, ehigh=90.0):
        '''cut_bins should be a list of N energy bin edges.
        cut_effs should be a list of N-1 cut efficiencies.
        '''
        if elow < 10.0:
            raise ValueError("Can't go lower than 10 MeV!")
        self.elow = elow
        self.ehigh = ehigh

        bins, effs = self._adjust_effs(cut_bins, cut_effs, elow, self.ehigh)
        self.cut_bins = bins
        self.cut_effs = effs
        self.efficiency_func = efficiency_func

        self.func = self._get_func(sknum)
        self.norm0 = self._get_norm0()
        self.norm = self._get_norm() # Normalization after cuts to 1

    def _adjust_effs(self, bins, effs, elow, ehigh):
        newbins, neweffs = bins, effs
        if elow < bins[0]:
            newbins, neweffs = [elow] + list(bins), [1.0] + list(effs)
        if bins[-1] < ehigh:
            newbins, neweffs = list(bins) + [ehigh], list(effs) + [1.0]
        return newbins, neweffs

    def _get_func(self, sknum):
        spall_emax = 24
        expcoeff = [3588.821, 1.60559875e5, 7768.694, 2181.575]
        exppow = [3.49809, 4.554, 3.7738, 3.3457]
        #exppow = [3.932, 4.650, 7.656, 3.746]
        #return lambda x: exp(-x * 0.836311) if x < spall_emax else 0
        return lambda x: exp(-(x**exppow[sknum - 1])/expcoeff[sknum - 1]) if x < spall_emax else 0

    def _get_norm0(self):
        ''' PDF normalization, before cuts '''
        spall_emax = 24
        return 1./quad(self.func, self.elow, min(self.ehigh,spall_emax))[0]

    def pdf_before_cuts(self, energy, region, ntag):
        ''' Properly normalized pdf, before cuts '''
        self._check_valid_energy(energy)
        if region != 1 or ntag: return 0
        return self.func(energy) * self.norm0

    def _get_norm(self):
        ''' Get pdf normalization after cuts '''
        area = 0.0
        for i, cut_eff in enumerate(self.cut_effs):
            e_lo, e_hi = self.cut_bins[i], self.cut_bins[i+1]
            if e_hi < self.elow:
                continue
            elif e_lo < self.elow:
                e_lo = self.elow
            if e_lo > self.ehigh:
                continue
            elif e_hi > self.ehigh:
                e_hi = self.ehigh
            a0, _ = quad(lambda x: self.pdf_before_cuts(x, 1, 0) * 
                         self.efficiency_func(x), e_lo, e_hi)
            area += a0 * cut_eff
        return 1.0 / area

    def _check_valid_energy(self, energy):
        if energy < self.elow or energy > self.ehigh:
            raise ValueError("Energy (%0.2f) outside range (%0.2f-%0.2f)"
                             % (energy, self.elow, self.ehigh))

    def pdf(self, energy, region, ntag = 0):
        ''' Properly normalized pdf, after cuts '''
        try:
            self._check_valid_energy(energy)
        except ValueError:
            return 1e-10

        ebin = digitize(energy, self.cut_bins)
        if ebin < 1 or ebin >= len(self.cut_bins):
            return 1e-10

        cut_eff = self.cut_effs[ebin-1]
        p0 = self.pdf_before_cuts(energy, region, ntag)
        return self.norm * cut_eff * p0 * self.efficiency_func(energy)


class relic_sk:
    '''
    Wrapper for relic signal pdf. For SK-IV, pdf has 3 Cherenkov angle
    regions and 2 ntag regions.
    '''
    def __init__(self, sknum, spectrum_energies, spectrum_values, ch_frac,
                 efficiency_func, efficiency_func_n=None, elow=16.0,
                 elow_n=None, ehigh=90.0, ntag_ebins=None, ntag_effs=None,
                 ntag_bgs=None, ntag_eff_ps=None, ntag_bg_ps=None):
        ''' efficiency_func should only depend on energy.'''
        if elow < 16.0:
            raise ValueError("Can't go lower than 16 MeV!")
        if len(spectrum_energies) != len(spectrum_values):
            raise ValueError("lengths of spectrum_energies and spectrum_values don't match")
        if elow < min(spectrum_energies):
            raise ValueError("elow must be within given spectrum")
        if elow_n is not None:
            if elow_n < min(spectrum_energies):
                raise ValueError("elow must be within given spectrum")

        self.sknum = sknum
        self.energies = spectrum_energies
        self.spec_values = spectrum_values
        self.effs = [efficiency_func, efficiency_func_n]
        self.elows = [elow, elow_n]
        self.ehigh = ehigh
        self.cherenkov_frac = ch_frac
        self.nregions_frac = [1.0]
        self.ntag_ebins = ntag_ebins
        self.ntag_effs = ntag_effs
        self.ntag_bgs = array(ntag_bgs)
        self.ntag_bg_ps = ntag_bg_ps

        if sknum == 4:
            ntagvars = [efficiency_func_n, ntag_ebins, ntag_effs,
                        ntag_bgs, ntag_eff_ps, ntag_bg_ps]
            for nv in ntagvars:
                assert nv is not None
            if elow_n is None:
                self.elows[1] = elow
            self.ntag_effs = array(ntag_effs) * ntag_eff_ps
            self.nregions_frac = self._get_ntag_fracs()

        self.spec_min = min(self.energies)
        self.spec = interpolate.interp1d(self.energies, self.spec_values,
                                         bounds_error=False, fill_value=0)
        self.norm0 = self._get_norm0() # Norm of source pdf to 1 in ana range
        self.norm_16_90 = self._get_norm_lims(16, 90) # Norm in 16-90 MeV range
        # self.norm_all = self._get_norm_all()  # Norm in whole range
        self.norm = self._get_norm() # Normalization after cuts to 1

    def _check_valid_energy(self, energy, ntag=False):
        if energy < self.elows[ntag] or energy > self.ehigh:
            raise ValueError("Energy (%0.2f) outside analysis range (%0.2f-%0.2f)"
                              % (energy, self.elows[ntag], self.ehigh))

    def _ntag_weight(self, E, ncap, ntag=True):
        ebin = digitize(E, self.ntag_ebins)
        ef = self.ntag_effs[ebin-1]
        bgrate = self.ntag_bgs[ebin-1]

        # Probability of tagging exactly 1 true neutron
        prob_1n = ncap * ef * (1 - ef)**(ncap - 1)
        # Probability of mistagging exactly 1 accidental
        prob_1b = self.ntag_bg_ps * bgrate * (1 - bgrate)**(self.ntag_bg_ps - 1)
        prob_0n = (1 - ef)**(ncap)
        prob_0b = (1 - bgrate)**(self.ntag_bg_ps)

        weight_1n = prob_1n * prob_0b + prob_1b * prob_0n
        if ntag:  # 1-neutron region
            return weight_1n
        else:  # 0 or multiple neutrons region
            return 1 - weight_1n

    def _get_ntag_fracs(self):
        fracs = [[], []]
        for en in self.ntag_ebins[:-1]:
            fracs[0] += [self._ntag_weight(en + 0.1, 1, ntag=False)]
            fracs[1] += [self._ntag_weight(en + 0.1, 1, ntag=True)]
        return fracs

    def _spectrum_before_cuts(self, energy, region, ntag=False):
        ''' Unnormed pdf '''
        ch_frac = self.cherenkov_frac[region]
        res = self.spec(energy) * ch_frac
        if self.sknum == 4:
            ebin = digitize(energy, self.ntag_ebins) - 1
            n_frac = self.nregions_frac[ntag][ebin]
            res *= n_frac
        return res

    def _get_norm0(self):
        area = 0.0
        for ntag in range(len(self.nregions_frac)):
            for region in range(3):
                area += quad(self._spectrum_before_cuts, self.elows[ntag],
                             self.ehigh, args=(region, ntag))[0]
        return 1.0 / area

    def _get_norm_lims(self, lo, hi):
        if lo < self.spec_min:
            raise ValueError("Outside interpolation range")
        area = 0.0
        for ntag in range(len(self.nregions_frac)):
            for region in range(3):
                area += quad(self._spectrum_before_cuts, lo, hi,
                             args=(region, ntag))[0]
        return 1.0 / area

    # def _get_norm_all(self):
    #     area = quad(self.spec, 0, 100)[0]
    #     return 1.0 / area

    def pdf_before_cuts(self, energy, region, ntag=False):
        ''' Properly normalized pdf, before cuts '''
        try:
            self._check_valid_energy(energy, ntag)
        except ValueError:
            return 1e-10

        return self._spectrum_before_cuts(energy, region, ntag) * self.norm0

    def _spectrum_after_cuts(self, energy, region, ntag=False):
        eff = self.effs[ntag](energy)
        # print(eff)
        ch_frac = self.cherenkov_frac[region]
        res = self.spec(energy) * eff * ch_frac
        if self.sknum == 4:
            ebin = digitize(energy, self.ntag_ebins) - 1
            n_frac = self.nregions_frac[ntag][ebin]
            res *= n_frac
        # print(energy, self.spec(energy))
        # print(res)
        return res

    def _get_norm(self):
        area = 0.0
        for ntag in range(len(self.nregions_frac)):
            for region in range(3):
                area += quad(self._spectrum_after_cuts, self.elows[ntag],
                             self.ehigh, args=(region, ntag))[0]
        return 1.0 / area

    def _get_ntag_sysfact(self):
        ''' Get multiplying factor for ntag systematics on overall efficiency '''
        area1 = 0
        area2 = 0
        for region in range(3):
            area2 += quad(self._spectrum_after_cuts, self.elows[1],
                         self.ehigh, args=(region, 1))[0]
        area1 = quad(lambda x: sum(self.spec(x) * self.effs[0](x)
                             * self.cherenkov_frac[r] 
                             * self.nregions_frac[1][digitize(x, self.ntag_ebins) - 1]
                             for r in range(3)), 
                             self.elows[0], self.ehigh)[0]
        return (area2- area1) * self.norm/self.norm0


    def overall_efficiency(self):
        ''' Cut efficiency within our energy region '''
        return self.norm0 / self.norm

    def overall_efficiency_16_90(self):
        ''' Cut efficiency relative to 16-90 MeV energy spectrum'''
        # print(self.norm_16_90)
        # print(self.norm)
        return self.norm_16_90 / self.norm

    # def overall_efficiency_all(self):
    #     ''' Cut efficiency over the whole srn spectrum '''
    #     return self.norm_all / self.norm

    def pdf(self, energy, region, ntag=False):
        ''' Properly normalized pdf, after cuts '''
        try:
            self._check_valid_energy(energy, ntag)
        except ValueError:
            return 1e-10

        return self.norm * self._spectrum_after_cuts(energy, region, ntag)
