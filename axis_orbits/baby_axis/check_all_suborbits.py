import sys
import os
import time
from collections import defaultdict, OrderedDict
from random import randint, sample
import numpy as np
from argparse import ArgumentParser
import shelve
from multiprocessing import Pool

from mmgroup import MM0, XLeech2, GCode, Octad, Xsp2_Co1
from mmgroup.general import Random_Subgroup

from mat22_orbits import SHELVE_NAME, load_orbits 
from mat22_orbits import BabyAxis 


sys.path.append(os.path.join("..", "utilities"))
from utilities import order_Nx0, is_Nx0_odd
from utilities import MM_to_GAP
from utilities_gap import create_input_for_gap, run_gap
sys.path.pop()


MP = 1   # Use multiprocessing for checks, when set

MAT22_SIZE = 22*21*20*16*3
N_X0_SIZE = (MAT22_SIZE * 2) << 34

OMEGA = XLeech2(0x800000)
NEG_OMEGA = -OMEGA

AXES = BabyAxis.representatives()


GAP_DIR = os.path.split(__file__)[0]
GAP_INPUT = os.path.join(GAP_DIR, "Nx0_orbit_structure.g")
GAP_OUTPUT = os.path.join(GAP_DIR, "Nx0_orbit_structure.txt")


def load_tables():
    global MAP_SUBORBIT, SUBORBIT_REPRESENTATIVES, SUBORBIT_SIZES
    global ORBITS, GENERATOS
    global SUBORBIT_CENTRALIZERS
    with shelve.open(SHELVE_NAME) as db:
        MAP_SUBORBIT = db["MAP_SUBORBIT"]
        SUBORBIT_REPRESENTATIVES = db["SUBORBIT_REPRESENTATIVES"]
        SUBORBIT_SIZES = db["SUBORBIT_SIZES"]
        SUBORBIT_CENTRALIZERS = db["SUBORBIT_CENTRALIZERS"]  
        ORBITS = load_orbits()



def suborbit_axis(axis):
    from watermark_suborbits import watermark_axis
    wm = watermark_axis(axis)
    return MAP_SUBORBIT[wm]



def test_hash_axes():
    from watermark_suborbits import suborbit_sample_axes
    hashes = set()
    for axis in suborbit_sample_axes():
        h_list = [axis.profile_Nxyz(t = (0,i))[1] for i in (0,1)]
        for h in h_list:
            assert h not in hashes
        for h in h_list:
            hashes.add(h)




def check_one_suborbit(i, ref_axis, entry, v, suborbit_size, centralizer, map):
    # orbit_name, entry, v = SUBORBIT_REPRESENTATIVES[i] 
    # ref_axis = BabyAxis.representatives()[orbit_name]
    # suborbit_size = int(SUBORBIT_SIZES[i])
    # centralizer = SUBORBIT_CENTRALIZERS[i]
    # map = MAP_SUBORBIT
    from watermark_suborbits import suborbit_sample_axes
    axis = suborbit_sample_axes()[i]
    assert axis == ref_axis *  Xsp2_Co1('c', v) ** -1
    e, f, s = order_Nx0(centralizer)
    order = f << (sum(e) + s)
    assert suborbit_size * order == N_X0_SIZE
    assert s == is_Nx0_odd(centralizer[0])
    assert axis * Xsp2_Co1(centralizer[0]) == axis
    for c in sample(centralizer[1:], 5):
        assert is_Nx0_odd(c) == 0
        assert axis * Xsp2_Co1(c) == axis
    for j in range(3):
        axis1 = axis * Xsp2_Co1('r', 'N_x0 & B')
        wm = (axis1.profile_Nxyz((0, 0), 1)[1],
               axis1.profile_Nxyz((0, 1), 1)[1])
        i1 = map[wm]
        assert i == i1, (i, i1)
    return e, f, s


def check_all_suborbits():
    list_cases = []
    map = MAP_SUBORBIT
    axis_representatives = BabyAxis.representatives()
    for i, (orbit_name, entry, v) in enumerate(SUBORBIT_REPRESENTATIVES):
        ref_axis = axis_representatives[orbit_name]
        suborbit_size = int(SUBORBIT_SIZES[i])
        centralizer = SUBORBIT_CENTRALIZERS[i]
        list_cases.append((
            i, ref_axis, entry, v, suborbit_size, centralizer, map))
    if MP:
        nprocesses = max(1, min(16, os.cpu_count() - 2))
        with Pool(processes = nprocesses) as pool:
            c = pool.starmap(check_one_suborbit, list_cases, chunksize = 8)
        pool.join()
    else:
        c = [check_one_suborbit(*args) for args in list_cases]
    with shelve.open(SHELVE_NAME) as db:
        db["SUBORBIT_SIZES_2"] = c




def display_suborbits():
    from watermark_suborbits import suborbit_sample_axes
    sample_axes = suborbit_sample_axes()
    hd_fmt = "%5s: %4s   %10s  %20s  %4s"
    print(hd_fmt % ("Orbit", "G_x0", "Images", "Orbit size", "#odd"))
    fmt = "%5d: %4s   {%3d, %3d}  %20d  %4s"
    with shelve.open(SHELVE_NAME) as db:
        suborbit_sizes2 = db["SUBORBIT_SIZES_2"]
    for i, _ in enumerate(SUBORBIT_REPRESENTATIVES):
        axis = sample_axes[i]
        orbit = axis.axis_type()
        assert suborbit_axis(axis) == i
        images = [suborbit_axis(axis * MM0('t', e)) for e in (1,2)]
        images.sort()
        size = SUBORBIT_SIZES[i]
        s = 2 >> suborbit_sizes2[i][2]
        print(fmt % (i, orbit, images[0], images[1], size, s))
    

def check_suborbits(check = True, verbose = False):
    load_tables()
    if check:
        #print(11)
        check_all_suborbits()
        #print(22)
        test_hash_axes()
    if verbose:
        display_suborbits()




def parse_args():
    description = ('Check centralizers of representaives of N_x0 '
    'orbits of 2A axes. ' 
    )
    usage = "usage: %prog [options]"
    parser = ArgumentParser(usage, description = description)
    parser.add_argument("-n",  dest="no_check", action="store_true",
        help = "Do not chckt orbits")
    parser.add_argument("-s",  dest="sizes", action="store_true",
        help = "Display sizes of orbits in Mat24")
    parser.add_argument("-g",  dest="gap", action="store_true",
        help = "Display input for GAP stucture description")
    options  = parser.parse_args()
    return options


def print_input_for_gap():
    f_in = GAP_INPUT
    f_out = GAP_OUTPUT
    with shelve.open(SHELVE_NAME) as db:
        centralizers = db["SUBORBIT_CENTRALIZERS"]
        orders = [x[1] for x in db["SUBORBIT_SIZES_2"]] 
    create_input_for_gap(centralizers, orders, f_in)

def call_gap():
    run_gap(GAP_INPUT, GAP_OUTPUT)


def ___display_output_for_gap():
    """Deprecated"""
    with shelve.open(SHELVE_NAME) as db:
        centralizers = db["SUBORBIT_CENTRALIZERS"]
    for i, c in enumerate(centralizers):
        g = MM_to_GAP(c)
        print(f"""G{i} = Group(
{g});
Print(G{i} = StructureDescription(G{i}, "\\n");
""")

if __name__ == "__main__":
    opt = parse_args()
    check_suborbits(check = not opt.no_check, verbose = opt.sizes)
    if opt.gap:
         display_output_for_gap()
