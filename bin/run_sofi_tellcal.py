import matplotlib as mpl
mpl.use('Agg')
mpl.rc("xtick", direction="in", labelsize=16)
mpl.rc("ytick", direction="in", labelsize=16)
mpl.rc("xtick.major", width=1., size=8)
mpl.rc("ytick.major", width=1., size=8)
mpl.rc("xtick.minor", width=1., size=5)
mpl.rc("ytick.minor", width=1., size=5)
import matplotlib.pyplot as plt

import os
import numpy as np
from astropy.io import fits
from glob import glob
import subprocess
from astropy.nddata import CCDData
from astropy.stats import sigma_clipped_stats
from sofitools import run_reduction
from sofitools import telluric_correction, find_extract_box, get_wavelength_arc
import argparse

parser = argparse.ArgumentParser(description='Calibrate the telluric feature')
parser.add_argument('-s', '--sci', required=True, dest='sci', type=str, 
                    help='Science 1D spectrum data')
parser.add_argument('-c', '--cal', required=True, dest='cal', type=str, 
                    help='Telluric standard data')
parser.add_argument('-t', '--type', default='A0V', dest='type', type=str, 
                    help='The stellar type of the telluric standard [A0V]')
parser.add_argument('-p', '--plot', dest='plot', default=False, action='store_true',
                    help='Plot the calibrated spectrum [False]')
parser.add_argument('--ncut', dest='ncut', default=20, type=int, 
                    help='Number of pixels to cut [20]')
args = parser.parse_args()


work_path = '.'
sci_file = '{0}/{1}'.format(work_path, args.sci)
cal_file = '{0}/{1}'.format(work_path, args.cal)
assert os.path.isfile(sci_file), 'Cannot find the science file ({})!'.format(args.sci)
assert os.path.isfile(cal_file), 'Cannot find the telluric file ({})!'.format(args.cal)

sci = fits.open(sci_file)
cal = fits.open(cal_file)
sci_filter = sci[0].header.get('ESO INS FILT1 NAME')
cal_filter = cal[0].header.get('ESO INS FILT1 NAME')
assert sci_filter == cal_filter, 'Science and telluric data are not in the same filter!'

wave = sci[1].data['Wavelength']
flux_sci = sci[1].data['Flux']
wave_cal = cal[1].data['Wavelength']
flux_cal = cal[1].data['Flux']

if (wave != wave_cal).any():
    flux_cal = np.interp(wave, wave_cal, flux_cal)

if args.type == 'A0V':
    tel_corr = telluric_correction.telluric_correction_A0V(wave, flux_cal)
elif args.type == 'B9V':
    tel_corr = telluric_correction.telluric_correction_B9V(wave, flux_cal)
else:
    raise Exception('The stellar model ({}) is not available!'.format(args.type))

flux_sci_corr = flux_sci * tel_corr

c1 = fits.Column(name='Wavelength', format='Float64', unit='micron', array=wave)
c2 = fits.Column(name='Flux', format='Float64', unit='adu', array=flux_sci_corr)
hdul = fits.HDUList([fits.PrimaryHDU(header=sci[0].header),
                     fits.BinTableHDU.from_columns([c1, c2])])
hdul[0].header['ESO PRO CATG'] = 'OBS_SPC1D_CORR'
calibrated_name = '{0}/{1}_calibrated.fits'.format(work_path, args.sci[:-5])
hdul.writeto(calibrated_name, overwrite=True)

if args.plot:
    fig_name = calibrated_name.split('/')[-1][:-5]
    obj_name = hdul[0].header.get('OBJECT')
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.step(wave, flux_sci_corr)
    ax.text(0.01, 0.05, '{0}, {1}'.format(obj_name, sci_filter), fontsize=24,
            transform=ax.transAxes, va='bottom', ha='left')
    
    fltr_wave = ((wave < 1.8) | (wave > 1.95)) & (wave > wave[args.ncut]) & (wave < wave[-args.ncut])
    mean, median, stddev = sigma_clipped_stats(flux_sci_corr[fltr_wave])
    ymin = median - 5 * stddev
    ymax = np.max(flux_sci_corr[fltr_wave] + 5 * stddev)
    
    ax.set_ylim([ymin, ymax])
    ax.set_xlabel('Wavelength ($\mu$m)', fontsize=24)
    ax.set_ylabel('Flux (ADU)', fontsize=24)
    ax.minorticks_on()
    plt.savefig('{0}/figs/{1}.pdf'.format(work_path, fig_name), bbox_inches='tight')
    plt.close()
