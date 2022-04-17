import numpy as np
from astropy.io import fits
from glob import glob
import subprocess
import os

cwd = os.getcwd()

def find_arc_raw(file_path):
    '''
    Find arc files.
    '''
    if not os.path.isabs(file_path):
        file_path = '{0}/{1}'.format(cwd, file_path)
    fList = sorted(glob('{}/*.fits'.format(file_path)))
    
    arcList = []
    for f in fList:
        header = fits.getheader(f, ext=0)
        dprcatg = header.get('ESO DPR CATG')
        dprtype = header.get('ESO DPR TYPE')
        if (dprcatg == 'CALIB') & (dprtype == 'LAMP'):
            arcList.append(f)
    return arcList
    

def find_flat_raw(file_path):
    '''
    Find flat files.
    '''
    if not os.path.isabs(file_path):
        file_path = '{0}/{1}'.format(cwd, file_path)
    fList = sorted(glob('{}/*.fits'.format(file_path)))
    
    flatList = []
    for f in fList:
        header = fits.getheader(f, ext=0)
        dprcatg = header.get('ESO DPR CATG')
        dprtype = header.get('ESO DPR TYPE')
        dprtech = header.get('ESO DPR TECH')
        if dprtech == 'IMAGE':
            continue
            
        if (dprcatg == 'CALIB') & (dprtype == 'FLAT'):
            flatList.append(f)
    return flatList
    

def find_science_raw(file_path):
    '''
    Find science files.
    '''
    if not os.path.isabs(file_path):
        file_path = '{0}/{1}'.format(cwd, file_path)
    fList = sorted(glob('{}/*.fits'.format(file_path)))
    
    sciDict = {}
    for f in fList:
        header = fits.getheader(f, ext=0)
        dprcatg = header.get('ESO DPR CATG')
        dprtype = header.get('ESO DPR TYPE')
        dprtech = header.get('ESO DPR TECH')
        if dprtech == 'IMAGE':
            continue
        
        name = header.get('ESO OBS TARG NAME')
        
        fltname =  header.get('ESO INS FILT1 NAME')
        assert fltname is not None, 'The filter is wrong??'
        if fltname not in sciDict:
            sciDict[fltname] = {}
        d = sciDict[fltname]
            
        if (dprcatg == 'SCIENCE') & (dprtype == 'OBJECT'):
            if name in d:
                d[name].append(f)
            else:
                d[name] = [f]
    return sciDict


def find_arc_reduced(work_path, filter, reduced_dir='reduced'):
    '''
    Find the reduced arc data.
    '''
    if not os.path.isabs(work_path):
        work_path = '{0}/{1}/{2}'.format(cwd, work_path, reduced_dir)
    else:
        work_path = '{0}/{1}'.format(work_path, reduced_dir)
    fList = sorted(glob('{}/*.fits'.format(work_path)))
    
    arcList = []
    for f in fList:
        header = fits.getheader(f, ext=0)
        procatg = header.get('ESO PRO CATG')
        fltname =  header.get('ESO INS FILT1 NAME')
        if (procatg == 'ARC_COEF') & (fltname == filter):
            arcList.append(f)
    return arcList
    

def find_flat_reduced(work_path, filter, reduced_dir='reduced'):
    '''
    Find the reduced flat data.
    '''
    if not os.path.isabs(work_path):
        work_path = '{0}/{1}/{2}'.format(cwd, work_path, reduced_dir)
    else:
        work_path = '{0}/{1}'.format(work_path, reduced_dir)
    fList = sorted(glob('{}/*.fits'.format(work_path)))
    
    flatList = []
    for f in fList:
        header = fits.getheader(f, ext=0)
        procatg = header.get('ESO PRO CATG')
        fltname =  header.get('ESO INS FILT1 NAME')
        if (procatg == 'MASTER_SP_FLAT') & (fltname == filter):
            flatList.append(f)
    return flatList


def find_science_combined(work_path, reduced_dir='reduced'):
    '''
    Find the reduced science data.
    '''
    if not os.path.isabs(work_path):
        work_path = '{0}/{1}/{2}'.format(cwd, work_path, reduced_dir)
    else:
        work_path = '{0}/{1}'.format(work_path, reduced_dir)
    fList = sorted(glob('{}/*.fits'.format(work_path)))
    
    sciDict = {}
    for f in fList:
        header = fits.getheader(f, ext=0)
        procatg = header.get('ESO PRO CATG')
        fltname =  header.get('ESO INS FILT1 NAME')
        assert fltname is not None, 'The filter is wrong??'
        
        if (procatg == 'OBS_COMBINED'):
            if fltname in sciDict:
                sciDict[fltname].append(f)
            else:
                sciDict[fltname] = [f]
    return sciDict
        

def write_arc_sof(file_list, work_path, calib_path, reduced_dir='reduced'):
    '''
    Write the SOF file of arc data.
    '''
    if not os.path.isabs(work_path):
        work_path = '{0}/{1}'.format(cwd, work_path)
    
    assert len(file_list) > 0, 'No data is found!'
    
    sof_path = '{0}/{1}/arc.sof'.format(work_path, reduced_dir)
    with open(sof_path, 'w') as f:
        for fn in file_list:
            f.write('{} SP_ARC\n'.format(fn))
        f.write('{}/xe.fits CALPRO_XE_CATALOG\n'.format(calib_path))
        f.write('{}/ne.fits CALPRO_NE_CATALOG\n'.format(calib_path))
    return sof_path
            

def write_flat_sof(file_list, work_path, reduced_dir='reduced'):
    '''
    Write the SOF file of flat data.
    '''
    if not os.path.isabs(work_path):
        work_path = '{0}/{1}'.format(cwd, work_path)
    
    assert len(file_list) > 0, 'No data is found!'
        
    sof_path = '{0}/{1}/flat.sof'.format(work_path, reduced_dir)
    with open(sof_path, 'w') as f:
        for fn in file_list:
            f.write('{} SP_FLAT\n'.format(fn))
    return sof_path
    
    
def write_jitter_sof(file_list, arc_list, flat_list, work_path, calib_path, 
                     reduced_dir='reduced'):
    '''
    Write the SOF file of science data.
    '''
    if not os.path.isabs(work_path):
        work_path = '{0}/{1}'.format(cwd, work_path)
    
    assert len(file_list) > 0, 'No data is found!'
        
    filename = file_list[0].split('/')[-1][:-5]
    sof_path = '{0}/{1}/{2}_jitter.sof'.format(work_path, reduced_dir, filename)

    with open(sof_path, 'w') as f:
        for fn in file_list:
            f.write('{} SP_NODDINGOBJ\n'.format(fn))
        for arc in arc_list:
            f.write('{} ARC_COEF\n'.format(arc))
        for flat in flat_list:
            f.write('{} MASTER_SP_FLAT\n'.format(flat))
        f.write('{}/oh.fits CALPRO_OH_CATALOG\n'.format(calib_path))
    return sof_path
