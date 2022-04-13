import numpy as np
from astropy.io import fits
from glob import glob
import subprocess

def find_arc(file_path):
    '''
    Find arc files.
    '''
    fList = sorted(glob('{}/*.fits'.format(file_path)))
    
    arcList = []
    for f in fList:
        header = fits.getheader(f, ext=0)
        dprcatg = header.get('ESO DPR CATG')
        dprtype = header.get('ESO DPR TYPE')
        if (dprcatg == 'CALIB') & (dprtype == 'LAMP'):
            arcList.append(f)
    return arcList
    

def find_flat(file_path):
    '''
    Find flat files.
    '''
    fList = sorted(glob('{}/*.fits'.format(file_path)))
    
    flatList = []
    for f in fList:
        header = fits.getheader(f, ext=0)
        dprcatg = header.get('ESO DPR CATG')
        dprtype = header.get('ESO DPR TYPE')
        if (dprcatg == 'CALIB') & (dprtype == 'FLAT'):
            flatList.append(f)
    return flatList
    

def find_science(file_path):
    '''
    Find science files.
    '''
    fList = sorted(glob('{}/*.fits'.format(file_path)))
    
    sciDict = {}
    for f in fList:
        header = fits.getheader(f, ext=0)
        dprcatg = header.get('ESO DPR CATG')
        dprtype = header.get('ESO DPR TYPE')
        name = header.get('ESO OBS TARG NAME')
        if (dprcatg == 'SCIENCE') & (dprtype == 'OBJECT'):
            if name in sciDict:
                sciDict[name].append(f)
            else:
                sciDict[name] = [f]
    return sciDict
    

def write_sof(namelist, datapath, datatype='sci'):
    '''
    '''
    
    

