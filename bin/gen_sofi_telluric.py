import os
import numpy as np
from sofitools import p2tool as p2t
import argparse

parser = argparse.ArgumentParser(description='Generate the OB of the telluric star')
parser.add_argument('-n', '--name', required=True, dest='name', type=str, 
                    help='Name of the telluric star')
parser.add_argument('-u', '--username', dest='username', type=str, 
                    help='ESO account user name')
parser.add_argument('-p', '--password', dest='password', type=str, 
                    help='ESO account password')
parser.add_argument('-r', '--runid', dest='runid', type=str, 
                    help='Observation run ID')
parser.add_argument('-f', '--folder', dest='folder', type=str, default='tmp',
                    help='Folder name')
args = parser.parse_args()

sofi = p2t.p2api_SOFI(args.runid, args.username, args.password, no_warning=True)
p2t.create_OB_telluric(args.name, sofi, args.folder)
