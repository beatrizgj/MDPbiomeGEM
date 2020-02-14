#!/usr/bin/python3

# Script to shape the desired output to be processed (MMODES)
# the datatable way

# @author: Jorge Carrasco Muriel
# Creation: 09/06/2019

import os
import re
import numpy as np
import datatable as dt
from datatable import f

def log(cons, media):
    '''
    Writes information of consortium object to file
    '''
    logf = 'simulations.txt'
    p = re.compile(r'#+ SIMULATION (\d+) #+')
    if os.path.isfile(logf): # parse last simulation number
        with open(logf) as l:
            for line in l.readlines():
                num_sim = p.search(line)
                if num_sim:
                    head = " SIMULATION "+str(int(num_sim.group(1))+1)+" "
    else:
        head = " SIMULATION 1 "
    lines = '{:{fill}{align}{width}}'.format(head,
        fill = '#',
        align = '^',
        width = 30) + "\n"
    lines += cons.__str__()
    pers = ', '.join([per["PERTURBATION"] for per in media])
    lines += "\nPERTURBATIONS: " + pers + "\n\n"
    with open(logf, "a") as l:
        l.write(lines)
    return

def equidistant(df, n):
    sample = np.linspace(df.nrows-1,1,n).astype('int')
    sample.sort()
    return df[sample, :]

def tsv_filter(medium = "", flux = "", txpers = {}, inplace = False, v = 0, equif = True, bin = False):
    '''
    Function that filters medium and fluxes TSV files based on perturbation times.
    INPUTS ->  medium: string, path to medium file;
            flux: string, path to medium file;
            txpers: dictionary, time : perturbation;
            inplace: bool, whether overwrite input paths (default False);
            v: float, volume magnitude to obtain medium concentrations;
            equif: bool, whether write an additional fluxes filtered file,
                    with 100 equidistant points (default True)
    OUTPUT -> it returns None, writes 2(3) TSV files
    '''
    dfs = []
    if not medium:
        print("Medium parameter wasn't supplied, it won't be generated.")
    else:
        dfs.append([dt.fread(medium), medium, 0])
        if v != 0:
            for i in range(1,dfs[0][0].ncols): dfs[0][0][:,i] = dfs[0][0][:,f[i]/v]
    if not flux:
        print("Medium parameter wasn't supplied, it won't be generated.")
    else:
        dfs.append([dt.fread(flux), flux, 1])
    if not medium:
        print("You must supply a txpers parameter. Exitting function...")
        return

    for log, path, n in dfs:
        log[:,'Perturbations'] = "FALSE" # now last column (-1)
        log[-1,-1] = "END"
        if len(txpers) > 1:
            for tp, per in txpers.items():
                if tp == 0:
                    log[0,-1] = per
                else:
                    # take last time that matches <= perturbation time
                    log[f.time == log[f.time < tp, f.time][-1,-1], -1] = per                
#                if per == 'START':
#                    log[0,-1] = 'START'
#                else:
#                    # take last index that matches <= perturbation time
#                    log[f.time == log[f.time <= tp, f.time][-1,-1], -1] = per
        else:
            log[0, -1] = 'START'
        if n != 0 and equif:
            log_equif = equidistant(log,100) # take 100 equidistant rows
            log_equif.to_csv(path[:-4] + '_equi' + '.tsv')
            del(log_equif)
            # TODO: I don't know how to implement a condroll with datatable
            # We aren't currentyly using it, anyway
        log = log[f.Perturbations != "FALSE", :]
        if inplace:
            log.to_csv(path)
        else:
            log.to_csv(path[:-4] + '_filtered' + '.tsv')
