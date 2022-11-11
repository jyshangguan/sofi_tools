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
calib_path = '{}/calib'.format(pipe_path)

flatList = run_reduction.find_flat_raw(work_path)
flatsof = run_reduction.write_flat_sof(flatList, work_path, reduced_dir)

esorex = '{}/bin/esorex'.format(pipe_path)
recipes = '--recipe-dir={}/lib/esopipes-plugins'.format(pipe_path)
output = '--output-dir={0}/{1}'.format(work_path, reduced_dir)
log = '--log-dir={0}/{1}'.format(work_path, reduced_dir)
sof_file = '{}/flat.sof'.format(reduced_dir)

process = subprocess.Popen([esorex, recipes, output, log, 'sofi_spc_flat', sof_file], 
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

os.system('mv {0}/*.paf {0}/{1}'.format(work_path, reduced_dir))
