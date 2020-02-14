#!/usr/bin/python3

# Script to take the output generated on atr2.
# TODO: this script should accept command line parameters for compressing, taking
# just fluxes/medium/plot, help display, etc.

# Author: Jorge Carrasco Muriel
# Date: 06/02/2019
# e-mail: jorge.cmuriel@alumnos.upm.es

import os
import sys
from shutil import copy as cp
import tarfile

def make_tarfile(output_filename, source_dir):
    ''' Tar.gz from directory '''
    with tarfile.open(output_filename, "w:gz") as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))

dest_dir = sys.argv[1]
if dest_dir.endswith('/'): # set a clean path to dest directory
    dest_dir = dest_dir + "results"
else:
    dest_dir = dest_dir + "/results"

dirs = [f for f in os.listdir(".") if os.path.isdir(f) and f.startswith("atr2")]
curr_path = os.getcwd()
if not os.path.exists(dest_dir):
    os.mkdir(dest_dir)
i = 1
n = len(dirs)
for d in dirs:
    sys.stdout.write(f"\rWorking on {d}... ({i}/{n})          ")
    sys.stdout.flush()
    i += 1
    if os.path.isfile(f"{d}/cons.p") and os.path.getsize(f"{d}/cons.p") > 0:
        cp(f"{d}/plot.tsv", f"{dest_dir}/medium_{d}.tsv")
        cp(f"{d}/filtered_fluxes.tsv", f"{dest_dir}/fluxes_{d}.tsv")
        cp(f"{d}/{d}_plot.png", f"{dest_dir}/")
    else:
        print(f"Simulation on directory {d} wasn't succesful")

# make a tar.gz to send it
print("Compressing results...", end = " ")
make_tarfile("results.tar.gz", dest_dir)
print("DONE!")
