import os
import numpy as np
__all__ = ['modulepath', 'bandpass_dict', 'cal_mag', 'flux_calibration']

# Obtain the current path
pathList = os.path.abspath(__file__).split("/")
# Create the path to the gData module
modulepath = "/".join(pathList[0:-1])

bp_2mass_j = np.loadtxt('{}/../data/bandpass/2MASS_J.dat'.format(modulepath))
bp_2mass_h = np.loadtxt('{}/../data/bandpass/2MASS_H.dat'.format(modulepath))
bp_2mass_ks = np.loadtxt('{}/../data/bandpass/2MASS_Ks.dat'.format(modulepath))
bandpass_dict = {
    '2MASS J': bp_2mass_j,
    '2MASS H': bp_2mass_h,
    '2MASS Ks': bp_2mass_ks,
}
f0_dict = {
    # http://www.adamgginsburg.com/filtersets.htm
    '2MASS J': 1594,  # Jy 
    '2MASS H': 1024,  # Jy
    '2MASS Ks': 666.7 # Jy
}

def cal_mag(wave, flux, band):
    '''
    Calculate the averaged magnitude of the spectrum in a band.
    '''
    assert (band in bandpass_dict), 'Cannot find the band ({}) in bandpass_dict!'
    assert (band in f0_dict), 'Cannot find the band ({}) in f0_dict!'
    
    bp = bandpass_dict[band]
    f0 = f0_dict[band]
    wave_bp = bp[:, 0]
    rsr_bp = bp[:, 1]
    
    flux_interp = np.interp(wave_bp, wave, flux)
    favg = np.trapz(flux_interp * rsr_bp, wave_bp) / np.trapz(rsr_bp, wave_bp)
    mag = -2.5 * np.log10(favg / f0)
    return mag

def flux_calibration(wave, flux, mag, band):
    '''
    '''
    assert (band in bandpass_dict), 'Cannot find the band ({}) in bandpass_dict!'
    assert (band in f0_dict), 'Cannot find the band ({}) in f0_dict!'
    
    bp = bandpass_dict[band]
    f0 = f0_dict[band]
    wave_bp = bp[:, 0]
    rsr_bp = bp[:, 1]
    
    flux_interp = np.interp(wave_bp, wave, flux)
    f_avg = np.trapz(flux_interp * rsr_bp, wave_bp) / np.trapz(rsr_bp, wave_bp)
    f_true = f0 * 10**(-0.4 * mag)
    
    r_cal = f_true / f_avg
    return r_cal
    
    