import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.modeling import models, fitting
from astropy.visualization import (PercentileInterval, LogStretch, SqrtStretch,
                                   ImageNormalize)
from astropy.stats import sigma_clipped_stats

__all__ = ['find_extract_box', 'get_wavelength_arc']


def get_wavelength_arc(filename, npix=1024):
    '''
    Calculate the wavelength according to the ARC_COEF data.
    
    Parameters
    ----------
    filename : str
        The filename of the ARC_COEF data.
    npix : int (default: 1024)
        The number of pixels in the dispersion dimension.
    
    Returns
    -------
    wave : 1d array
        Wavelength, units: Angstrom.
    '''
    wl_coef = fits.getdata(filename, ext=1)['WL_coefficients']
    pix = np.arange(npix)
    wave = np.polyval(wl_coef[::-1], pix+1)[::-1]
    return wave


def find_extract_box(data, nstd=5, ncut=20, std_init=10, plot=False, axs=None):
    '''
    Find the extracting box boundaries.
    
    Parameters
    ----------
    data : array_like
        2D spectrum data.
    nstd : float (default: 5)
        Number of pixels to expand on each side from the center for the 
        extracting box.  The box width is about 2*nstd (considering rounding to 
        integer).
    ncut : int (default: 20)
        The number of pixels to cut on the spatial direction of the data.
    std_init : float (default: 10)
        Initial guess of the width of the point spread function to fit the 
        center of the source, units: pixel.
    plot : bool (default: False)
        Plot the results if True.
    
    Returns
    -------
    x0, x1 : float
        The left and right boundaries of the extracting box, units: pixel.
    '''
    space_profile = np.median(data, axis=0)
    
    if isinstance(ncut, int):
        y = space_profile[slice(ncut, -ncut)]
        nstart = ncut
    elif isinstance(ncut, tuple): 
        y = space_profile[slice(ncut[0], ncut[1])]
        nstart = ncut[0]
    else:
        raise TypeError('The type of ncut is not correct!')
    
    x_p = np.argmax(y)
    
    g_init = models.Gaussian1D(amplitude=y[x_p], mean=x_p, stddev=std_init)
    fit_g = fitting.LevMarLSQFitter()
    
    x = np.arange(len(y))
    g = fit_g(g_init, x, y)
    
    mean = g.mean + nstart
    std = g.stddev
    x0, x1 = int(np.floor(mean - nstd * std)), int(np.ceil(mean + nstd * std))
    
    if plot:
        if axs is None:
            fig = plt.figure(figsize=(16, 8))
            ax0 = fig.add_axes([0.025, 0.05, 0.45, 0.9])
            ax1 = fig.add_axes([0.55, 0.55, 0.40, 0.40])
            ax2 = fig.add_axes([0.55, 0.05, 0.40, 0.40])
            axs = [ax0, ax1, ax2]
        else:
            ax0, ax1, ax2 = axs
        lw = 2
        
        # Plot the 2D spectrum
        ax = ax0
        _, vmed, vstd = sigma_clipped_stats(data, sigma=3)
        vmin = vmed - vstd
        vmax = vmed + 10 * vstd
        norm = ImageNormalize(data, vmin=vmin, vmax=vmax, stretch=LogStretch())
        ax.imshow(data, norm=norm, origin='lower')
        ax.axvspan(xmin=x0, xmax=x1, ls='--', fc='none', ec='C3', lw=lw)
        ax.set_xlabel('Space (pixel)', fontsize=24)
        ax.set_ylabel('Dispersion (pixel)', fontsize=24)
        ax.set_title('(a) 2D spectrum', fontsize=24)
        
        ax = ax1
        ax.plot(space_profile, label='Spatial profile')
        ax.plot(x + nstart, g(x), label='Fit')
        ax.axvspan(xmin=x0, xmax=x1, ls='--', fc='none', ec='C3', lw=lw, label='Extraction')
        ax.legend(loc='best', fontsize=16)
        ax.set_ylabel('Intensity (ADU)', fontsize=24)
        ax.set_title('(b) Spatial collapsed', fontsize=24)
        ax.minorticks_on()
        
        ax = ax2
        ax.plot(y, label='Spatial profile')
        ax.plot(x, g(x), label='Fit')
        ax.axvspan(xmin=x0-nstart, xmax=x1-nstart, ls='--', fc='none', ec='C3', lw=lw, label='Extraction')
        ax.set_xlabel('Space (pixel)', fontsize=24)
        ax.set_ylabel('Intensity (ADU)', fontsize=24)
        ax.set_title('(c) Spatial zoomed to fit', fontsize=24)
        ax.minorticks_on()
        return x0, x1, axs
        
    return x0, x1