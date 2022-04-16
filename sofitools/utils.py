import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy import units as u
from astroplan import FixedTarget
from astroplan import Observer
from astroplan.plots import plot_airmass 

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
    
    
def cal_airmass(ra_target, dec_target, ra_telluric, dec_telluric, time, observer, delta_time, plot=False, ax=None):
    '''
    Calculate the airmasses.
    '''
    c_targ = SkyCoord(ra=ra_target, dec=dec_target, unit=(u.hourangle, u.deg))
    c_tell = SkyCoord(ra=ra_telluric, dec=dec_telluric, unit=(u.hourangle, u.deg))
    
    t_targ = FixedTarget(coord=c_targ, name=None)
    airmass_target = observer.altaz(time, t_targ).secz.value
    
    t_tell = FixedTarget(coord=c_tell, name=None)
    time_tell = time + delta_time * u.hour
    airmass_tell = observer.altaz(time_tell, t_tell).secz.value
    
    if plot:
        observe_time = time + np.linspace(0, 15, 55)*u.hour
        
        if ax is None:
            fig, ax = plt.subplots(figsize=(15, 8))
        plot_airmass(t_targ, observer, observe_time, ax=ax, brightness_shading=True, style_kwargs=dict(color='C3'))
        plot_airmass(t_tell, observer, observe_time, ax=ax, brightness_shading=True, style_kwargs=dict(color='C0'))
        
        ax.axvline(x=time.plot_date, ls='--', color='C3')
        ax.axhline(y=airmass_target, ls='--', color='C3')
        ax.axvline(x=time_tell.plot_date, ls='--', color='C0')
        ax.axhline(y=airmass_tell, ls='--', color='C0')
        ax.plot([time_tell.plot_date, time_tell.plot_date], [airmass_target, airmass_tell], ls='-', lw=5, color='yellow')
        t = '\n'.join(['Target: {0}, {1}'.format(ra_target, dec_target), 
                       'Telluric: {0}, {1}'.format(ra_telluric, dec_telluric)])
        ax.text(0.01, 1.01, t, fontsize=18, transform=ax.transAxes, va='bottom', ha='left')
        t = '\n'.join(['Separation: {0:.1f} degree'.format(c_targ.separation(c_tell).degree),
                       'delta_airmass: {0:.1e}'.format(np.abs(airmass_target-airmass_tell))])
        ax.text(0.99, 1.01, t, fontsize=18, transform=ax.transAxes, va='bottom', ha='right')
        t = '{0} (+{1} hour)'.format(time.iso, delta_time)
        ax.text(0.02, 0.04, t, fontsize=18, transform=ax.transAxes, va='bottom', ha='left')
        ax.minorticks_on()
    return airmass_target, airmass_tell