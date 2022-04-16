import os
import numpy as np
from astropy.io import fits
from glob import glob
import subprocess
from sofitools import run_reduction

work_path = './'
reduced_dir = 'reduced'
if not os.path.isdir(reduced_dir):
    os.makedirs(reduced_dir)
    
pipe_path = os.environ['SOFIPIPELINE']
calib_path = '{}/calib/sofi-1.5.13'.format(pipe_path)
sciDict = run_reduction.find_science_raw(work_path)

esorex = '{}/bin/esorex'.format(pipe_path)
recipes = '--recipe-dir={}/lib/esopipes-plugins'.format(pipe_path)
output = '--output-dir={0}/{1}'.format(work_path, reduced_dir)
log = '--log-dir={0}/{1}'.format(work_path, reduced_dir)

if len(sciDict) == 0:
    raise RuntimeError('Cannot find any data!')
else:
    print('#----------------------------')
    print('# Find data in these filters: {}'.format(list(sciDict.keys())))
    print('#----------------------------\n')

for k_fltr in sciDict:
    arcList = run_reduction.find_arc_reduced(work_path, k_fltr, reduced_dir)
    flatList = run_reduction.find_flat_reduced(work_path, k_fltr, reduced_dir)
    for tn in sciDict[k_fltr]:
        sciList = sciDict[k_fltr][tn]
        scisof = run_reduction.write_jitter_sof(sciList, arcList, flatList, work_path, 
                                                calib_path, reduced_dir)
        combined_name = '{0}_combined.fits'.format(scisof[:-4])
        
        if os.path.isfile(combined_name):
            print('{}\n  File exists, skip the reduction...\n'.format(combined_name))
            continue
        
        process = subprocess.Popen([esorex, recipes, log, 'sofi_spc_jitter', scisof], 
                                   stdout=subprocess.PIPE,
                                   universal_newlines=True)
        
        while True:
            output = process.stdout.readline()
            print(output.strip())
            # Do something else
            return_code = process.poll()
            if return_code is not None:
                print('RETURN CODE', return_code)
                # Process has finished, read rest of the output 
                for output in process.stdout.readlines():
                    print(output.strip())
                break
        
        process = subprocess.Popen(['mv', '{0}/sofi_spc_jitter_combined.fits'.format(work_path), combined_name],
                     stdout=subprocess.PIPE, 
                     stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()

os.system('mv {0}/*.paf {0}/{1}'.format(work_path, reduced_dir))
