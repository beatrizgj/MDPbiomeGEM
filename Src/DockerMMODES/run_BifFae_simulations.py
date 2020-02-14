#!/usr/bin/python3

# Script to generate permutations of media and run a MMODES simulation in different
# parallel processes.

import sys, os, re
import random
import json
import mmodes
import dill as pickle
import subprocess
import pandas as pd
from copy import deepcopy as dcp
from multiprocessing import Pool, Lock
sys.path.append("/home/jorge/mmgit/")
from data_gen import *

lock = Lock() # global variable, accessed by child processes

def runn(md, all_random = False):
    global lock
    lock.acquire()
    print("Starting simulation on", md)
    lock.release()

    subprocess.run(["mkdir", md])
    subprocess.run(["cp", "-r", "fb1/ModelsInput", md])
    path_curr = os.getcwd()
    os.chdir(md)
    mod_file = "ModelsInput"
    intervl = 8

    files = [f for f in os.listdir(".") if f.find("log_template.txt") != -1 or f[-3:] == "tsv"]
    subprocess.run(["rm"]+files)
    with open(mod_file+"/media.json") as json_file:
        gen_media = json.load(json_file)

    inulin = {
        "PERTURBATION" : "Inulin",
        "MEDIA" : {
            "Inulin_29FruGlc" : "0.0005"
        }
    }
    fos = {
        "PERTURBATION" : "FOS",
        "MEDIA" : {
            "Kestose_C18H32O16" : "0.005"
        }
    }
    starch = {
        "PERTURBATION" : "Starch",
        "MEDIA" : {
            "Starch_11Glc" : "0.00128"
        }
    }
    faepraa2165 = {
        "PERTURBATION" : "Pro-Fprausx3.5",
        "MEDIA" : {
            "FAEPRAA2165" : "0.00035"
        }
    }

    media = [gen_media[0]]
    pers_perm = [inulin, fos, starch]
    if all_random:
        for i in range(0, 5):
            if random.random() > 0.75:
                pertadd = dcp(faepraa2165)
                rand = [0.00001, 0.0001]
            else:
                pertadd = dcp(random.choice(pers_perm))
                if pertadd == fos:
                    rand = [0.005,0.015]
                elif pertadd == inulin:
                    rand = [0.0002,0.0006]
                else:
                    rand = [0.0005, 0.0015]
            perval = random.uniform(rand[0], rand[1])
            pertadd["MEDIA"][list(pertadd["MEDIA"])[0]] = perval
            pertadd["PERTURBATION"] += str(round(perval,5))
            media.append(pertadd)
    else:
        for i in range(0, 5):
            if random.random() > 0.9:
                media.append(faepraa2165)
            else:
                media.append(random.choice(pers_perm))

    fp = random.uniform(0.0001, 0.00055)
    ba = random.uniform(0.0001, 0.00055)
    # 1) instantiate Consortium
    volume_petri = 2*3.141593*45*45*15*1e-6 # 2pi*r²*h -> mm³ to L
    cons = mmodes.Consortium(stcut = -10, mets_to_plot = ["Kestose_C18H32O16", "Inulin_29FruGlc", "Starch_11Glc"],
    v = volume_petri, manifest = "COMETS_manifest.txt", work_based_on = "name", max_growth = 10, title=md)
    # 2) add models
    cons.add_model(mod_file+"/iFap484.V01.00.xml", float(fp), solver = "glpk", method = "fba")
    cons.add_model(mod_file+"/iBif452.V01.00.xml", float(ba), solver = "glpk", method = "fba")# dMets = {glc.id: glc})

    lock.acquire()
    log(cons, media)
    lock.release()

    it = 1
    t_pers = []
    pers = []
    for mper in media:
        if it == 1:
            # 4) instantiate media
            cons.media = cons.set_media(mper["MEDIA"], True)
        else:
            # 4.2) or add perturbation
            for k in mper["MEDIA"]:
                if mper["MEDIA"][k] == 0:
                    cons.media[k] = 0
            cons.add_mets(mper["MEDIA"], True)
        pers.append(mper["PERTURBATION"])
        # 5) run it
        t_pers.append(cons.T[-1])
        cons.run(verbose = False, plot = False, maxT = intervl+cons.T[-1],
        integrator = "vode", stepChoiceLevel = (0.,0.5,100000.), outp = md+"_plot.png")
        it += 1

    txpers = {t_pers[i]: pers[i] for i in range(len(t_pers))}
    # 6. Save relevant objects and generate output files
    pickle.dump(cons, open("cons.p", "wb"))
    pickle.dump(txpers, open("txpers.p", "wb"))
    pickle.dump(media, open("media.p", "wb"))
    #mmodes.vis.plot_comm(cons)
    gen_biomass("plot.tsv", cons.v, cons.T[:-1], txpers)
    gen_fluxes("flux_log_template.txt", cons.T, cons.manifest.models)
    os.chdir(path_curr)

    lock.acquire()
    print(f"\033[1;32;40mProcess with directory {md} out!\033[0m")
    del(cons)
    del(txpers)
    del(media)
    lock.release()

    return

def main():
    if len(sys.argv) > 1:
        try:
            cl = int(sys.argv[1])
            battery = ["fb2","fb3", "fb4", "fb5", "fb6"]
        except:
            print("You need to provide one integer as first argument >= 1 or none.")
            sys.exit(1)
        if len(sys.argv) > 2:
            battery = ["fb"+str(i+2) for i in range(int(sys.argv[2]))]
    else:
        print("Running on 1 process as default...")
        runn("debugfb"); sys.exit(0)
    # work in FCFS queue
    piscine = Pool(processes = cl)
    piscine.map(runn, battery)
    # work in batch
    # i = 0
    # while i < len(battery):
    #     piscine = Pool(processes = cl) # workaround accumulation of memory in workers
    #     for j in range(cl):
    #         if i >= len(battery):
    #             break
    #         piscine.apply_async(runn, args=(battery[i],))
    #         i += 1
    #     piscine.close()
    #     piscine.join()

if __name__ == "__main__":
    main()
