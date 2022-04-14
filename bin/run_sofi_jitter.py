import os
import numpy as np
from astropy.io import fits
from glob import glob
import subprocess
from sofitools import run_reduction

work_path = './'
pipe_path = os.environ['SOFIPIPELINE']
calib_path = '{}/calib/sofi-1.5.13'.format(pipe_path)
esorex = '{}/bin/esorex'.format(pipe_path)

sciDict = run_reduction.find_science_raw(work_path)

if len(sciDict) == 0:
    raise RuntimeError('Cannot find any data!')
else:
    print('#----------------------------')
    print('# Find data in these filters: {}'.format(sciDict.keys()))

for k_fltr in sciDict:
    arcList = run_reduction.find_arc_reduced(work_path, k_fltr)
    flatList = run_reduction.find_flat_reduced(work_path, k_fltr)
    
    for tn in sciDict[k_fltr]:
        sciList = sciDict[k_fltr][tn]
        scisof = run_reduction.write_jitter_sof(sciList, arcList, flatList, work_path=work_path, calib_path=calib_path)
        
        process = subprocess.Popen([esorex, 'sofi_spc_jitter', scisof], 
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
        
        process = subprocess.Popen(['mv', '{}/sofi_spc_jitter_combined.fits'.format(work_path), 
                                    '{0}_combined.fits'.format(scisof[:-4])],
                     stdout=subprocess.PIPE, 
                     stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
