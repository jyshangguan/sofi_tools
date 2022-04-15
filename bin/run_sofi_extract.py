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
from sofitools import run_reduction
from sofitools import find_extract_box, get_wavelength_arc
import argparse

parser = argparse.ArgumentParser(description='Extract 1D spectrum')
parser.add_argument('-w', '--nstd', default=5, dest='nstd', type=float, 
                    help='Extraction box (half) width')
parser.add_argument('-c', '--ncut', default=20, dest='ncut', type=int, 
                    help='Cut the edges of the detector')
parser.add_argument('-i', '--std_init', default=10, dest='std_init', type=float, 
                    help='Initial guess of the Gaussian fitting')
parser.add_argument('-p', '--plot', default=False, dest='plot', 
                    action='store_true', help='Plot the extraction details if True')
args = parser.parse_args()


work_path = '.'
pipe_path = os.environ['SOFIPIPELINE']
calib_path = '{}/calib/sofi-1.5.13'.format(pipe_path)
esorex = '{}/bin/esorex'.format(pipe_path)
recipes = '--recipe-dir={}/lib/esopipes-plugins'.format(pipe_path)

sciDict = run_reduction.find_science_combined(work_path)

if len(sciDict) == 0:
    raise RuntimeError('Cannot find any data!')
else:
    print('#----------------------------')
    print('# Find data in these filters: {}'.format(list(sciDict.keys())))
    print('#----------------------------\n')

for k_fltr in sciDict:
    arcList = run_reduction.find_arc_reduced(work_path, k_fltr)
    wave = get_wavelength_arc(arcList[0]) / 1e4  # Micron
    
    sciList = sciDict[k_fltr]
    for sci_name in sciList:
        spc1d_name = '{}_spc1d.fits'.format(sci_name[:-5])
        
        if os.path.isfile(spc1d_name):
            print('{}\n  File exists, skip the reduction...\n'.format(spc1d_name))
            continue
        
        sci = CCDData.read(sci_name, unit='adu')
        
        if args.plot:
            if not os.path.exists('{}/figs'.format(work_path)):
                os.makedirs('{}/figs'.format(work_path))
                
            fig = plt.figure(figsize=(14, 14))
            ax0 = fig.add_axes([0.02, 0.40, 0.45, 0.45])
            ax1 = fig.add_axes([0.55, 0.65, 0.40, 0.2])
            ax2 = fig.add_axes([0.55, 0.40, 0.40, 0.2])
            ax3 = fig.add_axes([0.02, 0.04, 0.93, 0.3])
            axs = (ax0, ax1, ax2)
        else:
            axs = None
                
        res = find_extract_box(sci, nstd=args.nstd, ncut=args.ncut, 
                               std_init=args.std_init, plot=args.plot, 
                               axs=axs)
        sci_spc1d = sci.data[:, res[0]:res[1]].sum(axis=1)
        
        c1 = fits.Column(name='Wavelength', format='Float64', unit='micron', array=wave)
        c2 = fits.Column(name='Flux', format='Float64', unit='adu', array=sci_spc1d)
        hdul = fits.HDUList([fits.PrimaryHDU(header=sci.header),
                             fits.BinTableHDU.from_columns([c1, c2])])
        hdul[0].header['ESO PRO CATG'] = 'OBS_SPC1D'
        hdul.writeto(spc1d_name, overwrite=True)
        
        if args.plot:
            fig_name = spc1d_name.split('/')[-1][:-5]
            obj_name = hdul[0].header.get('OBJECT')
            
            ax3.step(wave, sci_spc1d)
            ax3.text(0.01, 0.05, '{0}, {1}'.format(obj_name, k_fltr), fontsize=24,
                     transform=ax3.transAxes, va='bottom', ha='left')
            ax3.set_xlabel('Wavelength ($\mu$m)', fontsize=24)
            ax3.set_ylabel('Flux (ADU)', fontsize=24)
            ax3.minorticks_on()
            plt.savefig('{0}/figs/{1}.pdf'.format(work_path, fig_name), bbox_inches='tight')
            plt.close()
                        
                        
