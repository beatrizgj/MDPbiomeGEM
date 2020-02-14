#!/usr/bin/python3

# Script PHASE 2 of simulation with atrazine consortium.
# 1.  SIMULATION of 5h, with 4 possible perturbations, in  particular, adding  phosphate to the medium, or incrementing the biomass of H.stevensii and/or Halobacillus sp.
# Command: 'python3 run_Atrazine_simulations.py num_processors num_simulations'

# Author: Jorge Carrasco Muriel
# Date: 30/04/2019
# Modified: 04/11/2019 by Beatriz Garcia Jimenez

import matplotlib as MPL
MPL.use('Agg') # change backend of display. I'm not really sure if this will work on a container
import sys, os, re
import subprocess
import json
import mmodes
import random
import dill as pickle
from time import sleep
from multiprocessing import Process, Lock
from data_gen import tsv_filter, log
from copy import deepcopy

lock = Lock() # global variable, accessed by child processes

def runn(md = ""):
    lock.acquire()
    print("Starting simulation on", md)
    lock.release()

    # make working directory based on argument
    subprocess.run(["mkdir", md])
    # clean directory if needed
    files = [f'{md}/{f}' for f in os.listdir(".") if f.find("log_template.txt") != -1 or f[-3:] == "tsv" or f[-3:] == 'png']
    if files:
        subprocess.run(["rm"]+files)

    # PARAMETERS
    # 1) Common parameters of simulation
    media_file = "media.json"
    mod_dir = "/home/cleanapp/draft/ModelsInput/"
    intervl = 1 # time of simulation between perturbations; total time will be 2h
    
    # 2) random parameters of simulation
    # 2.1) Medium
    with open(mod_dir + media_file) as json_file:
        gen_media = json.load(json_file)
    # 2.2) Biomasses
    brand = lambda: random.uniform(0.000013, 0.000055)
    biomasses = [brand()] # BGJ: Biomass for athrobacter that always will be in the consortium
    # we want to take all the possibilities with *equal* probabilities
    # BGJ: add 2 additional values to the biomasses vector, that could be only one, both ones or any.
    chosen = random.random()
    if chosen > 0.75:
        biomasses += [brand(), 0]
    elif chosen > 0.5:
        biomasses += [0, brand()]
    elif chosen > 0.25:
        biomasses += [brand(), brand()]
    else:
        biomasses += [0,0]
    ar, hb, hl = biomasses # BGJ: split biomasses in 3 different variables
 
    # SIMULATION
    # 1) instantiate Consortium
    volume_petri = 2*3.141593*45*45*15*1e-6 # 2pi*r²*h -> mm³ to L
    cons = mmodes.Consortium(stcut = 1e-8, v = volume_petri, comets_output = False,
    manifest = f"{md}/fluxes.tsv", work_based_on = "id", max_growth = 10,
    mets_to_plot = ["cpd03959_e0", "cpd00027_e0"], title = "Atrazine "+md)

    # 2.1) add models, with random biomass
    cons.add_model(mod_dir+"Arthrobacter_CORRECTED.json", float(ar), solver = "glpk", method = "fba")
    cons.add_model(mod_dir+"Halobacillus_sp_CORRECTED.json", float(hb), solver = "glpk", method = "fba")
    cons.add_model(mod_dir+"Halomonas_stevensii_CORRECTED.json", float(hl), solver = "glpk", method = "fba")

    st_ids = [id for id in cons.models if not id.startswith("k")]
    print(f"Models were loaded on {md}.")

    # 2.2) define initial medium => medium or exudate, already in JSON file
    root = deepcopy(gen_media[1])
    del(gen_media[random.randint(0,1)]) # choose medium or exudate

    # 2.3) add several random perturbations : strain | both | none
    fix_biomass=0.000034  # BGJ: fix biomass to add in the perturbation as probiotics
    for pert in range(3): # BGJ: 2019.11.08
        # perturbations always carry atrazine
        gen_media.append({"MEDIA" : {"cpd03959_e0" : 0.15}, "PERTURBATION" : "ATRAZINE"})
        chosen = random.random()
        if chosen > 0.8:
            gen_media[-1]["MEDIA"][st_ids[0]] = fix_biomass*100
            gen_media[-1]["PERTURBATION"] = st_ids[0]
        elif chosen > 0.6: # 2019.11.12: BGJ: to increase biomass 2nd strain: st_ids[1]
            gen_media[-1]["MEDIA"][st_ids[1]] = fix_biomass*100
            gen_media[-1]["PERTURBATION"] = st_ids[1]
        elif chosen > 0.4:
            gen_media[-1]["PERTURBATION"] = st_ids[0] + "_" + st_ids[1]
            gen_media[-1]["MEDIA"][st_ids[0]] = fix_biomass*100
            gen_media[-1]["MEDIA"][st_ids[1]] = fix_biomass*100
   #    elif chosen > 0.2:
   #        gen_media[-1]["MEDIA"] = root["MEDIA"]
   #        gen_media[-1]["PERTURBATION"] = "ROOT_EXUDATE"
        else:
            gen_media[-1]["MEDIA"]['cpd00009_e0'] = 1
            gen_media[-1]["PERTURBATION"] = 'PHOSPHATE'       
       #else: none perturbation (only atrazine)

    # 2.3) add 2nd perturbation (nothing)
    # we need this to evaluate the last state of the consortium because the
    # reward function in MDPbiome for this particular case evaluates degradation of
    # atrazine 1 h after the simulation
    gen_media.append({"MEDIA" : {"cpd03959_e0" : 0.15}, "PERTURBATION" : "ATRAZINE"})

    # 3) Write log
    lock.acquire()
    log(cons, gen_media)
    lock.release()

    it = 1
    t_pers = []
    pers = []
    for mper in gen_media:
        if it == 1:
            # 4) instantiate media
            cons.media = cons.set_media(mper["MEDIA"], True)
        else:
            for k in mper["MEDIA"]: # not needed at all...
                if mper["MEDIA"][k] == 0:
                    cons.media[k] = 0
            cons.add_mets(mper["MEDIA"], True)
        pers.append(mper["PERTURBATION"])
        # 5) run it
        t_pers.append(cons.T[-1])
        cons.run(verbose=False, plot=False, maxT = intervl+cons.T[-1], integrator = "FEA",
        stepChoiceLevel=(0.00027,0.5,100000.), outp = f'{md}/{md}_plot.png', outf = f'{md}/plot.tsv')
        it += 1
    txpers = { t_pers[i]: pers[i] for i in range(len(t_pers)) }

    # 6. Save simulation
    with open(f'{md}/cons.p', 'wb') as f:
        pickle.dump(cons, f)
    with open(f'{md}/txpers.p', 'wb') as f:
        pickle.dump(txpers, f)
    tsv_filter(f'{md}/plot.tsv', f'{md}/fluxes.tsv', txpers, inplace= False, v = cons.v)
    if os.path.isfile(f'{md}/fluxes_filtered.tsv'):
        os.unlink(f'{md}/fluxes.tsv')
    mmodes.vis.plot_comm(cons)

    lock.acquire()
    print(f"\033[1;32;40mProcess with directory {md} out!\033[0m")
    lock.release()

    del(cons)
    del(txpers)
    del(gen_media)

    return

def main():
    if len(sys.argv) > 1:
        try:
            cl = int(sys.argv[1])
            num_sims = cl
        except:
            print("You need to provide one integer as first argument >= 1 or none.")
            sys.exit(1)
        if len(sys.argv) > 2:
            num_sims = int(sys.argv[2])
    else:
        print("Running on 1 process as default...")
        runn("debugrfb"); sys.exit(0)
    # don't overwrite previous simulations
    prev_sims = len([d for d in os.listdir(".") if d.startswith("atr2.")]) + 1
    battery = ["atr2."+str(prev_sims+i) for i in range(num_sims)]
    print(f'{len(battery)} tasks have been loaded in queue with {cl} cores.\n')
    my_pool_loop(runn, battery, cl)
    return

def my_pool_loop(f, all_args, cl = 6):
    '''
    Emulate pool with different processes + check if child is stuck
    INPUTS:
        - f: function to be parallelized
        - all_args: arguments to be consumed by the function
    '''
    dict_map = {} # pid : arguments, since they get inaccesible when child is started
    def rawcount(filename):
        '''
        Fast line counting pure Python
        https://stackoverflow.com/a/27518377/11272105
        '''
        with open(filename, 'rb') as f:
            lines = 0
            buf_size = 1024 * 1024
            read_f = f.raw.read
            buf = read_f(buf_size)
            while buf:
                lines += buf.count(b'\n')
                buf = read_f(buf_size)
        return lines

    def is_stuck(child):
        ''' Check whether output of a process has changed '''
        snap = rawcount(f'{dict_map[str(child.pid)]}/plot.tsv')
        sleep(10)
        snap_t = rawcount(f'{dict_map[str(child.pid)]}/plot.tsv')
        if snap != snap_t:
            return False
        else:
            return True

    battery = deepcopy(all_args)
    piscine = []
    i = 0
    while i < len(battery)-1:
        for _ in range(cl-len(piscine)):
            piscine.append(Process(target=f, args = (battery[i],)))
            piscine[-1].start()
            dict_map[str(piscine[-1].pid)] = battery[i]
            print(f'Process {piscine[-1].pid} spawned! (with args = {battery[i]})')
            i += 1
            added = True
            if i >= len(battery)-1:
                break
        if added:
            # This is important to output in this precise moment all the prints
            sys.stdout.write('Parent taking a nap just 1 min...\n')
            sys.stdout.flush()
            sleep(60) # wait for children to initialize and start their simulations
            added = False
        for j in range(len(piscine)):
            if not piscine[j].is_alive():
                piscine[j] = None
            elif is_stuck(piscine[j]):
                print(f'Process {piscine[j].pid} with args = {dict_map[str(piscine[j].pid)]} got stuck! Removing...')
                battery.append(dict_map[str(piscine[j].pid)]) # append argument to the queue
                piscine[j].terminate()
                piscine[j] = None
        piscine = [p for p in piscine if p != None]
    for p in piscine: # wait last tasks to end
        print('Last tasks have spawned! Going to sleep...')
        p.join()
    return

if __name__ == "__main__":
    main()
