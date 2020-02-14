#!/usr/bin/python3

# Script to take the output generated on atr2.

import os
import sys
from shutil import copy as cp
import tarfile
import mmodes
from mmodes.vis import plot_comm
import dill as pickle
from matplotlib import rcParams, rcParamsDefault, get_backend, rcParamsOrig

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
    prefix=d.split(".")
    #comm=pickle.load(open(f"{d}/cons.p", "rb"))
    #comm.outplot="{d}/{d}_plot.png"
    #plot_comm(comm)

    cp(f"{d}/fluxes_equi.tsv", f"{dest_dir}/equi_fluxes_{prefix[0]}_{prefix[1]}.tsv")
    cp(f"{d}/plot_filtered.tsv", f"{dest_dir}/biomass_{prefix[0]}_{prefix[1]}.tsv")
    cp(f"{d}/fluxes_filtered.tsv", f"{dest_dir}/fluxes_{prefix[0]}_{prefix[1]}.tsv")
    if os.path.isfile(f"{d}/{d}_plot.png"):
        cp(f"{d}/{d}_plot.png", f"{dest_dir}/{prefix[0]}_{prefix[1]}_plot.png")
    
# make a tar.gz to send it
make_tarfile("results.tar.gz", dest_dir)
print("DONE")
