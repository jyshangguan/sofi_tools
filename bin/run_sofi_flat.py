import os
import numpy as np
from astropy.io import fits
from glob import glob
import subprocess
from sofitools import run_reduction
from pathlib import Path

work_path = './'
pipe_path = os.environ['SOFIPIPELINE']
calib_path = '{}/calib/sofi-1.5.13'.format(pipe_path)
esorex = '{}/bin/esorex'.format(pipe_path)

flatList = run_reduction.find_flat_raw(work_path)
flatsof = run_reduction.write_flat_sof(flatList, work_path=work_path)

process = subprocess.Popen([esorex, 'sofi_spc_flat', 'flat.sof'], 
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

