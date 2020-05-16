from numpy import digitize
from pickle import load
from scipy import interpolate
from scipy.integrate import quad

class SK4pdf:
    ''' Usage: 
    nc = SK4pdf(2, cut_edges, efficiencies)
    nc.pdf(energy, region)
    '''
    def __init__(self, ev_type, cut_bins, cut_effs, pdf_dir):
        ''' 0 < ev_type < 3 for nue, numu, nc, or mupi background.
        cut_bins should be a list of N energy bin edges.
        cut_effs should be a list of N-1 cut efficiencies.
        '''
        self.ev_type = ev_type
        self.cut_bins = cut_bins
        self.cut_effs = cut_effs

        # Get background pdf spline fit parameters
        if ev_type == 0:  # nue
            self.tck = [self._read_tck('%s/cc_e%d.p'% (pdf_dir, i)) for i in range(3)]
        elif ev_type==1:  # numu
            self.tck = [self._read_tck('%s/cc_mu%d.p'% (pdf_dir, i)) for i in range(3)]
        elif ev_type==2:  # nc
            self.tck = [self._read_tck('%s/nc%d.p'% (pdf_dir, i)) for i in range(3)]
        elif ev_type==3:  # mupi
            self.tck = [self._read_tck('%s/mupi%d.p'% (pdf_dir, i)) for i in range(3)]
        self.norm0 = self._get_norm0() # Normalization of source pdf to 1
        self.norm = self._get_norm() # Normalization after cuts to 1

    def _read_tck(self, file):
        with open(file, 'r') as fl: tck = load(fl)
        return tck    
    
    def _pdf_before_cuts_raw(self, energy, region):
        ''' Unnormed pdf '''
        p = interpolate.splev(energy, self.tck[region], der=0)
        return p.clip(min=0) # Enforce non-negativity
    
    def _get_norm0(self):
        area = 0.0
        for region in range(3):
            area += quad(self._pdf_before_cuts_raw, 16, 90, args=(region,))[0]
        return 1.0 / area

    def pdf_before_cuts(self, energy, region):
        ''' Properly normalized pdf, before cuts '''
        return self._pdf_before_cuts_raw(energy, region) * self.norm0

    def _get_norm(self):
        area = 0.0
        for region in range(3):
            for i, cut_eff in enumerate(self.cut_effs):
                e_lo, e_hi = self.cut_bins[i], self.cut_bins[i+1]
                a0, _ = quad(self.pdf_before_cuts, e_lo, e_hi, args=(region,))
                area += a0 * cut_eff
        return 1.0 / area
    
    def pdf(self, energy, region):
        ''' Properly normalized pdf, after cuts '''
        ebin = digitize(energy, self.cut_bins)
        if ebin < 1 or ebin >= len(self.cut_bins):
            cut_eff = 1.0 # outside of cut bins
        else: cut_eff = self.cut_effs[ebin-1]
        p0 = self.pdf_before_cuts(energy, region)
        return self.norm * cut_eff * p0

