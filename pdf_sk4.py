'''
Wrapper classes for relic and SK-IV backgrounds.
'''
from pickle import load as loadpick
from numpy import digitize
from scipy import interpolate
from scipy.integrate import quad

class bg_sk4:
    ''' Usage:
    nc = bg_sk4(2, cut_edges, efficiencies)
    nc.pdf(energy, region)
    '''
    def __init__(self, ev_type, cut_bins, cut_effs, cut_bins_n, cut_effs_n,
                 pdf_dir, elow, elow_n=None):
        ''' 0 < ev_type < 3 for nue, numu, nc, or mupi background.
        cut_bins should be a list of N energy bin edges.
        cut_effs should be a list of N-1 cut efficiencies.
        '''
        if elow < 16.0:
            raise ValueError("Can't go lower than 16 MeV!") #TODO: change?
        if elow_n is None:
            elow_n = elow
        self.elows = [elow, elow_n]
        self.ehigh = 90.0

        bins, effs = self._adjust_effs(cut_bins, cut_effs, elow, self.ehigh)
        bins_n, effs_n = self._adjust_effs(cut_bins_n, cut_effs_n, elow_n, self.ehigh)
        self.cut_bins = [bins, bins_n]
        self.cut_effs = [effs, effs_n]
        self.ev_type = ev_type
        self.pdf_dir = pdf_dir

        self.params = self._load_param_array(self.ev_type)
        self.norm0 = self._get_norm0() # Normalization of source pdf to 1
        self.norm = self._get_norm() # Normalization after cuts to 1


    def _adjust_effs(self, bins, effs, elow, ehigh):
        newbins, neweffs = bins, effs
        if elow < bins[0]:
            newbins, neweffs = [elow] + bins, [1.0] + effs
        if bins[-1] < ehigh:
            newbins, neweffs = bins + [ehigh], effs + [1.0]
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
            # return 1e-10

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
                    a0, _ = quad(self.pdf_before_cuts, e_lo,
                                 e_hi, args=(region, ntag))
                    area += a0 * cut_eff
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
        return self.norm * cut_eff * p0

class relic_sk:
    '''
    Wrapper for relic signal pdf. For SK-IV, pdf has 3 Cherenkov angle
    regions and 2 ntag regions.
    '''
    def __init__(self, sknum, spectrum_energies, spectrum_values,
                 efficiency_func, efficiency_func_n=None,
                 elow=16.0, elow_n=None):
        ''' efficiency_func should only depend on energy.'''
        if elow < 16.0:
            raise ValueError("Can't go lower than 16 MeV!")
        if len(spectrum_energies) != len(spectrum_values):
            raise ValueError("lengths of spectrum_energies and spectrum_values don't match")
        if elow < min(spectrum_energies):
            raise ValueError("elow must be within given spectrum")

        self.sknum = sknum
        self.energies = spectrum_energies
        self.spec_values = spectrum_values
        # self.eff = efficiency_func
        # self.eff_n = efficiency_func_n
        self.effs = [efficiency_func, efficiency_func_n]
        self.elows = [elow, elow_n]
        self.ehigh = 90.0
        self.cherenkov_frac = [9.433e-04, 9.925e-01, 6.525e-03]
        self.nregions_frac = [1.0]

        if sknum == 4:
            assert efficiency_func_n is not None
            if elow_n is None:
                self.elows[1] = elow
            self.nregions_frac = [0.1, 0.9] # TODO: Update neutron fractions?

        self.spec = interpolate.interp1d(self.energies, self.spec_values,
                                         bounds_error=False, fill_value=0)
        self.norm0 = self._get_norm0() # Normalization of source pdf to 1
        self.norm = self._get_norm() # Normalization after cuts to 1

    def _check_valid_energy(self, energy, ntag=False):
        if energy < self.elows[ntag] or energy > self.ehigh:
            raise ValueError("Energy (%0.2f) outside range (%0.2f-%0.2f)"
                              % (energy, self.elows[ntag], self.ehigh))

    def _spectrum_before_cuts(self, energy, region, ntag=False):
        ''' Unnormed pdf '''
        n_frac = self.nregions_frac[ntag]
        ch_frac = self.cherenkov_frac[region]
        return self.spec(energy) * ch_frac * n_frac

    def _get_norm0(self):
        area = 0.0
        for ntag in range(len(self.nregions_frac)):
            for region in range(3):
                area += quad(self._spectrum_before_cuts, self.elows[ntag],
                             self.ehigh, args=(region, ntag))[0]
        return 1.0 / area

    def pdf_before_cuts(self, energy, region, ntag=False):
        ''' Properly normalized pdf, before cuts '''
        try:
            self._check_valid_energy(energy, ntag)
        except ValueError:
            return 1e-10

        return self._spectrum_before_cuts(energy, region, ntag) * self.norm0

    def _spectrum_after_cuts(self, energy, region, ntag=False):
        eff = self.effs[ntag](energy)
        n_frac = self.nregions_frac[ntag]
        ch_frac = self.cherenkov_frac[region]
        return self.spec(energy) * eff * n_frac * ch_frac

    def _get_norm(self):
        area = 0.0
        for ntag in range(len(self.nregions_frac)):
            for region in range(3):
                area += quad(self._spectrum_after_cuts, self.elows[ntag],
                             self.ehigh, args=(region, ntag))[0]
        return 1.0 / area

    def overall_efficiency(self):
        ''' Cut efficiency over the entire pdf '''
        return self._get_norm0() / self._get_norm()

    def pdf(self, energy, region, ntag=False):
        ''' Properly normalized pdf, after cuts '''
        try:
            self._check_valid_energy(energy, ntag)
        except ValueError:
            return 1e-10

        return self.norm * self._spectrum_after_cuts(energy, region, ntag)
