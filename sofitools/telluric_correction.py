import numpy as np
from .utils import modulepath

__all__ = ['spec_a0v', 'spec_vega', 'telluric_correction_A0V', 'telluric_correction_B9V']


def telluric_correction_A0V(wave, flux):
    '''
    Calculate the telluric correction spectrum for an A0V standard star.
    This is still a very preliminary method.
    
    Parameters
    ----------
    wave, flux : 1d array
        Wavelength and flux of the standard star.
    
    Returns
    -------
    tel_corr : 1d array
        Telluric correction spectrum.
    '''
    telspec = spec_a0v
    wave_tel = telspec['wave']
    flux_tel = telspec['flux']
    
    flux_tel_interp = np.interp(wave, wave_tel, flux_tel)
    tel_corr = flux_tel_interp / flux
    tel_corr /= np.median(tel_corr)
    return tel_corr
    

def telluric_correction_B9V(wave, flux):
    '''
    Calculate the telluric correction spectrum for an A0V standard star.
    This is still a very preliminary method.
    
    Parameters
    ----------
    wave, flux : 1d array
        Wavelength and flux of the standard star.
    
    Returns
    -------
    tel_corr : 1d array
        Telluric correction spectrum.
    '''
    telspec = spec_b9v
    wave_tel = telspec['wave']
    flux_tel = telspec['flux']
    
    flux_tel_interp = np.interp(wave, wave_tel, flux_tel)
    tel_corr = flux_tel_interp / flux
    tel_corr /= np.median(tel_corr)
    return tel_corr
    

templatepath = '{}/../data/telluric_template'.format(modulepath)

# From http://cdsarc.u-strasbg.fr/viz-bin/vizExec/Vgraph?J/PASP/110/863/./uka0v&
spec_a0v_pickles1998 = np.loadtxt('{}/Pickles1998_A0V.dat'.format(templatepath))
spec_a0v_pickles1998[:, 0] /= 1e3  # Convert to micron

fltr = (spec_a0v_pickles1998[:, 0] > 1.2) & (spec_a0v_pickles1998[:, 0] < 2.5)
wave = spec_a0v_pickles1998[fltr, 0]
flux = spec_a0v_pickles1998[fltr, 1]

contWR = [(1.38, 1.43), (1.69, 1.72), (1.75, 1.80), (2.01, 2.15), (2.19, 2.4)]
fltr_cont = np.zeros_like(wave, dtype=bool)
for (w1, w2) in contWR:
    fltr = (wave>= w1) & (wave <= w2)
    fltr_cont[fltr] = True
wave_cont = wave[fltr_cont]
flux_cont = flux[fltr_cont]

p_cont = np.polyfit(wave_cont, flux_cont, deg=4)
flux_cont_m = np.polyval(p_cont, wave)
flux_norm = flux / flux_cont_m

spec_a0v = {
    'wave': wave,
    'flux': flux,
    'flux_cont_m': flux_cont_m,
    'flux_norm': flux_norm
}

spec_b9v_pickles1998 = np.loadtxt('{}/Pickles1998_B9V.dat'.format(templatepath))
spec_b9v_pickles1998[:, 0] /= 1e3  # Convert to micron
fltr = (spec_b9v_pickles1998[:, 0] > 1.2) & (spec_b9v_pickles1998[:, 0] < 2.5)
wave = spec_b9v_pickles1998[fltr, 0]
flux = spec_b9v_pickles1998[fltr, 1]
spec_b9v = {
    'wave': wave,
    'flux': flux,
}


# From: http://kurucz.harvard.edu/stars/vega/; I manually cut only the NIR part.
spec_vega_dat = np.loadtxt('{}/vegallpr25.50000resam5_cut_nir.dat'.format(templatepath))
spec_vega = {
    'wave': spec_vega_dat[:, 0],
    'flux': spec_vega_dat[:, 1],
}
