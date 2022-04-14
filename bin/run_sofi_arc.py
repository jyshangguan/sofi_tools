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
recipes = '--recipe-dir={}/lib/esopipes-plugins'.format(pipe_path)

arcList = run_reduction.find_arc_raw(work_path)
arcsof = run_reduction.write_arc_sof(arcList, work_path=work_path, calib_path=calib_path)

process = subprocess.Popen([esorex, recipes, 'sofi_spc_arc', 'arc.sof'], 
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

