import sys
import os
import time
from collections import defaultdict, OrderedDict
from random import randint, sample
import numpy as np
from argparse import ArgumentParser
import shelve

from mmgroup import MM0, MM, AutPL, Xsp2_Co1, XLeech2

from mat22_orbits import SHELVE_NAME 
from mat22_orbits import load_orbits, load_samples
from mat22_orbits import BabyAxis, configure_axis_group 
sys.path.append(os.path.join("..", "utilities"))
from utilities import trim_N_x0
sys.path.pop()

OMEGA = XLeech2(0x800000)
NEG_OMEGA = -OMEGA
CENTRALIZER_SIZE = 30

configure_axis_group()

DICT_NAME = "monster_tables"

AXES = BabyAxis.representatives()
ORBITS = list(AXES.keys())
#print(ORBITS)
ORBIT_KEYS = {}  # for sorting orbits
for orbit in ORBITS:
    ORBIT_KEYS[orbit] = (AXES[orbit].stage, int(orbit[:-2]), 
                           orbit[-2], 1 - int(orbit[-1]))

orbits_dict = load_orbits()
samples_dict = load_samples()
reps_dict, lengths_dict = {}, {}
for name, orbits in orbits_dict.items():
    reps, lengths = orbits.representatives()
    reps_dict[name] = reps
    lengths_dict[name] = lengths



def watermark_axis(axis):
    _, h_even, _ = axis.profile_Nxyz((0,0), 1)
    _, h_odd, _ = axis.profile_Nxyz((0,1), 1)
    return min(h_even, h_odd),  max(h_even, h_odd)


def sample_from_numpy_ints(a, n_samples):
    indices = sample(range(len(a)), min(len(a), n_samples))
    return [int(a[i]) for i in indices]

def suborbit_represetatives(axis_type):
    axis0 = AXES[axis_type]
    orbit_key = ORBIT_KEYS[axis_type]
    reps, lengths = orbits_dict[axis_type].representatives()
    samples = samples_dict[axis_type]
    assert len(reps) == len(samples) == len(lengths)
    for v, sample, length in zip(reps, samples, lengths):
        g_transform  = Xsp2_Co1('c', v) ** -1
        axis = axis0 * g_transform
        watermark0 = watermark_axis(axis)
        #print(hex(v), length)
        for d in [v1 for v1 in sample if v1]:
                axis_d = axis0 * MM0('c', d) ** -1 #  * MM0('l', 1)
                watermark_d = watermark_axis(axis_d) 
                assert watermark_d == watermark0
                axis_d1 = axis_d * MM0('r', 'N_x0 & B')
                watermark_d1 = watermark_axis(axis_d1) 
                assert watermark_d1 == watermark0
        # Compute centralizer of axis in N_x0
        c = [orbits_dict[axis_type].rand_stabilizer(v) ** g_transform
            for i in range(CENTRALIZER_SIZE)]
        c = trim_N_x0(c)
        axis_key =  orbit_key + (length,)
        #print(axis_key)
        yield axis, axis_key, c    




def check_watermarks():
    global WATERMARK_DICT
    S_OK = "All %2d suborbits could be separated"
    S_BAD = "Suborbit clusters of sizes %s found"
    bad_cases = []
    total_orbits = 0
    total_watermarks = set()
    watermark_dict = defaultdict(list)
    clusters = []
    n = 0
    for orbit in ORBITS:
        for i, (axis, key, c) in enumerate(suborbit_represetatives(orbit)):
            watermark = watermark_axis(axis)
            n += 1
            watermark_dict[watermark].append((orbit, i, axis, key, c))
            #print(watermark)
            total_watermarks.add(watermark) 
            for g in c:
                OMEGA_transformed = OMEGA * g 
                assert OMEGA_transformed in [OMEGA, NEG_OMEGA]
    ok = n == len(total_watermarks) == len(watermark_dict)
    if ok:
        print(S_OK % n)
        watermark_dict_new = {} 
        for w, value_list in watermark_dict.items():
            assert len(value_list) == 1
            watermark_dict_new[w] = value_list[0]
        return watermark_dict_new
    else:
        for w, value_list in watermark_dict.items():
            if len(value_list) > 1:
                clusters.append(len(value_list))
                bad_cases.append(value_list)
                print([(x[0], x[1], x[3]) for x in value_list]) 
        print(S_BAD % (clusters))
        ERR = "Watermarking of suborbits has failed"
        raise  ValueError(ERR)  

def suborbit_key(watermark):
    orbit, i, axis, axis_key, _ = WATERMARK_DICT[watermark]
    w_axes = [watermark_axis(axis * MM0('t', e)) for e in (1,2)]
    key1, key2 = [WATERMARK_DICT[w][3] for w in w_axes]
    if key1 > key2:
        key1, key2 = key2, key1
    return axis_key, key1, key2, watermark


def number_suborbits(verbose = True):
    global WATERMARK_DICT
    global MAP_SUBORBIT
    WATERMARK_DICT = check_watermarks()
    MAP_SUBORBIT = {}
    suborbits = list(WATERMARK_DICT.keys())
    suborbits = sorted(suborbits, key = suborbit_key)
    #print(len(suborbits))
    last = (None,)
    for i, watermark in enumerate(suborbits):
        current = suborbit_key(watermark)
        assert current != last
        if verbose:
            print(i, current[:-1])
            if current[:-1] == last[:-1]:
                print("Non-canonic disambiguation!")
        MAP_SUBORBIT[watermark] = i
        last = current

def enhance_map_suborbit():
    global MAP_SUBORBIT
    map_new = {}
    for (h0, h1), suborbit in MAP_SUBORBIT.items(): 
        map_new[(h1,h0)] = suborbit
    MAP_SUBORBIT.update(map_new)


def suborbit_to_representative(verbose = 1):
    global SUBORBIT_REPRESENTATIVES, SUBORBIT_SIZES
    global SUBORBIT_CENTRALIZERS
    N = 93150  // 2  #  No of type-4 vector in Leech lattice mod 2
    a = [None] * len(MAP_SUBORBIT)
    SUBORBIT_SIZES = [None] * len(MAP_SUBORBIT)
    SUBORBIT_CENTRALIZERS = [None] * len(MAP_SUBORBIT)
    with shelve.open(SHELVE_NAME) as db:
        orbit_sizes = db["ORBIT_SIZES"]
    for name in ORBITS:
        #reps, lengths = orbits_dict[name].representatives()
        reps, lengths = reps_dict[name], lengths_dict[name]
        axis0 = AXES[name]
        for i, data in enumerate(reps):
            axis = axis0 * MM0('c', data) ** -1
            watermark = watermark_axis(axis)
            suborbit_no = MAP_SUBORBIT[watermark]
            a[suborbit_no] = (name, i, data)
            orbit_size = orbit_sizes[name]
            size, mod = divmod(lengths[i] * orbit_size, N)
            assert mod == 0
            SUBORBIT_SIZES[suborbit_no] = size
            centralizer = WATERMARK_DICT[watermark][4] 
            SUBORBIT_CENTRALIZERS[suborbit_no] = [
                str(MM(g)) for g in centralizer ]
    for entry in a:
        assert entry is not None
    SUBORBIT_REPRESENTATIVES = a
    assert sum(SUBORBIT_SIZES) == sum(orbit_sizes.values())
    if verbose:
        for i, entry in enumerate(a):
            print(i, entry) 


def suborbit_sample_axes():
    with shelve.open(SHELVE_NAME) as db:
        representative_data = db["SUBORBIT_REPRESENTATIVES"]
    representatives = []
    for i, (orbit_name, entry, v) in enumerate(representative_data):
        assert v == reps_dict[orbit_name][entry]
        representatives.append(AXES[orbit_name] * MM0('c', v) ** -1) 
    return representatives


def final_check():
    from check_all_suborbits import check_suborbits
    print("Checking N_x0 of suborbits ...")
    check_suborbits(check = True, verbose = False)           
    print("Check of suborbits passed")

       
       
VERBOSE = 0   
def watermark_suborbits():
    number_suborbits(verbose = VERBOSE)
    suborbit_to_representative(verbose = VERBOSE)
    enhance_map_suborbit()
    with shelve.open(SHELVE_NAME) as db:
        db["MAP_SUBORBIT"] = MAP_SUBORBIT
        db["SUBORBIT_REPRESENTATIVES"] = SUBORBIT_REPRESENTATIVES
        db["SUBORBIT_SIZES"] = SUBORBIT_SIZES
        db["SUBORBIT_CENTRALIZERS"] = SUBORBIT_CENTRALIZERS
    print("Suborbit data written to shelve")
    final_check()



if __name__ == "__main__":
   watermark_suborbits()


