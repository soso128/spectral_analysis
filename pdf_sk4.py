from numpy import digitize
from pickle import load as loadpick
from scipy import interpolate
from scipy.integrate import quad

class bg_sk4:
    ''' Usage: 
    nc = bg_sk4(2, cut_edges, efficiencies)
    nc.pdf(energy, region)
    '''
    def __init__(self, ev_type, cut_bins, cut_effs, pdf_dir, elow, ntag=False):
        ''' 0 < ev_type < 3 for nue, numu, nc, or mupi background.
        cut_bins should be a list of N energy bin edges.
        cut_effs should be a list of N-1 cut efficiencies.
        '''
        if elow < 16.0: raise ValueError("Can't go lower than 16 MeV!")
        self.ev_type = ev_type
        self.cut_bins = cut_bins
        self.cut_effs = cut_effs
        self.elow = elow
        self.ehigh = 90.0
        if cut_bins[0] > self.elow: 
            self.cut_bins, self.cut_effs = [self.elow] + self.cut_bins, [1.0] + self.cut_effs
        if cut_bins[-1] < self.ehigh:
            self.cut_bins, self.cut_effs = self.cut_bins + [self.elow], self.cut_effs + [1.0]

        # Get background pdf spline fit parameters
        suf = '_ntag' if ntag else ''
        if ev_type == 0:  # nue
            self.tck = [self._read_tck('%s/cc_e%d%s.p'% (pdf_dir, i, suf)) for i in range(3)]
        elif ev_type==1:  # numu
            self.tck = [self._read_tck('%s/cc_mu%d%s.p'% (pdf_dir, i, suf)) for i in range(3)]
        elif ev_type==2:  # nc
            self.tck = [self._read_tck('%s/nc%d%s.p'% (pdf_dir, i, suf)) for i in range(3)]
        elif ev_type==3:  # mupi
            self.tck = [self._read_tck('%s/mupi%d%s.p'% (pdf_dir, i, suf)) for i in range(3)]
        self.norm0 = self._get_norm0() # Normalization of source pdf to 1
        self.norm = self._get_norm() # Normalization after cuts to 1

    def _read_tck(self, file):
        with open(file, 'rb') as fl: tck = loadpick(fl, encoding = "bytes")
        return tck    
    
    def _pdf_before_cuts_raw(self, energy, region):
        ''' Unnormed pdf '''
        p = interpolate.splev(energy, self.tck[region], der=0)
        return p.clip(min=0) # Enforce non-negativity
    
    def _get_norm0(self):
        area = 0.0
        for region in range(3):
            area += quad(self._pdf_before_cuts_raw, self.elow, self.ehigh, args=(region,))[0]
        return 1.0 / area
    
    def _check_valid_energy(self, energy):
        if energy < self.elow or energy > self.ehigh:
            #raise ValueError("Energy (%0.2f) outside range (%0.2f-%0.2f)" % (energy, self.elow, self.ehigh))
            return 1e-10
        return

    def pdf_before_cuts(self, energy, region):
        ''' Properly normalized pdf, before cuts '''
        self._check_valid_energy(energy)
        return self._pdf_before_cuts_raw(energy, region) * self.norm0

    def _get_norm(self):
        area = 0.0
        for region in range(3):
            for i, cut_eff in enumerate(self.cut_effs):
                e_lo, e_hi = self.cut_bins[i], self.cut_bins[i+1]
                if e_hi < self.elow: continue
                elif e_lo < self.elow: e_lo = self.elow
                a0, _ = quad(self.pdf_before_cuts, e_lo, e_hi, args=(region,))
                area += a0 * cut_eff
        return 1.0 / area
    
    def pdf(self, energy, region):
        ''' Properly normalized pdf, after cuts '''
        self._check_valid_energy(energy)
        ebin = digitize(energy, self.cut_bins)
        if ebin < 1 or ebin >= len(self.cut_bins): 
            #raise ValueError("This shouldn't happen...")
            return 1e-10
        cut_eff = self.cut_effs[ebin-1]
        p0 = self.pdf_before_cuts(energy, region)
        return self.norm * cut_eff * p0

class relic_sk4:
    def __init__(self, spectrum_energies, spectrum_values, efficiency_func, elow=16.0):
        ''' efficiency_func should only depend on energy.'''
        if elow < 16.0: raise ValueError("Can't go lower than 16 MeV!")
        if len(spectrum_energies) != len(spectrum_values): 
            raise ValueError("lengths of spectrum_energies and spectrum_values don't match")
        if elow < min(spectrum_energies):
            raise ValueError("elow must be within given spectrum")
        self.energies = spectrum_energies
        self.spec_values = spectrum_values
        self.eff = efficiency_func
        self.elow = elow
        self.ehigh = 90.0
        self.cherenkov_frac = [9.433e-04, 9.925e-01, 6.525e-03]

        self.spec = interpolate.interp1d(self.energies, self.spec_values, bounds_error = False, fill_value = 0)
        self.norm0 = self._get_norm0() # Normalization of source pdf to 1
        self.norm = self._get_norm() # Normalization after cuts to 1
    
    def _check_valid_energy(self, energy):
        if energy < self.elow or energy > self.ehigh:
            #raise ValueError("Energy (%0.2f) outside range (%0.2f-%0.2f)" % (energy, self.elow, self.ehigh))
            return 1e-10
        return

    def _spectrum_before_cuts(self, energy, region):
        ''' Unnormed pdf '''
        return self.spec(energy) * self.cherenkov_frac[region]
    
    def _get_norm0(self):
        area = 0.0
        for region in range(3):
            area += quad(self._spectrum_before_cuts, self.elow, self.ehigh, args=(region,))[0]
        return 1.0 / area
    
    def pdf_before_cuts(self, energy, region):
        ''' Properly normalized pdf, before cuts '''
        self._check_valid_energy(energy)
        return self._spectrum_before_cuts(energy, region) * self.norm0

    def _spectrum_after_cuts(self, energy, region):
        return self.spec(energy) * self.eff(energy) * self.cherenkov_frac[region]

    def _get_norm(self):
        area = 0.0
        for region in range(3):
            area += quad(self._spectrum_after_cuts, self.elow, self.ehigh, args=(region,))[0]
        return 1.0 / area

    def overall_efficiency(self):
        return self._get_norm0()/self._get_norm()
    
    def pdf(self, energy, region):
        ''' Properly normalized pdf, after cuts '''
        self._check_valid_energy(energy)
        return self.norm * self._spectrum_after_cuts(energy, region)
