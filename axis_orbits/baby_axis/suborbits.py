import sys
import os
import time
from collections import defaultdict, OrderedDict
from random import randint, sample
import numpy as np
from argparse import ArgumentParser
import shelve
from multiprocessing import Pool

from mmgroup import MM0, XLeech2, GCode, Octad
from mmgroup.general import Orbit_Lin2

from mat22_orbits import AXES, configure_axis_group
from mat22_orbits import SHELVE_NAME 
from mat22_orbits import load_orbits, load_samples

configure_axis_group() 



def centralizers(orbit_name):
    return orbits[orbit_name].generators()



def map_leech2_vector(v):
    return MM0('c', v) ** -1



T1 = MM0('t', 1); T2 = MM0('t', 2) 
def py_triality_orbits(axis, g = MM0()):
    ax = axis * g
    lst = [(ax * T1).axis_type(), (ax * T2).axis_type()]
    lst.sort()
    return tuple(lst)

triality_orbits = py_triality_orbits

def is_good_axis(axis):
    orbit = axis.axis_type()
    g = axis.reduce_G_x0()
    ref_axis = AXES[orbit]
    return (axis * g).v15 == ref_axis.v15



N = 8
N_VERIFY = 1

def process_orbit(orbit_name, axis, orbits, samples, verbose = False):
    d = defaultdict(int)
    if verbose:
        print("orbit", orbit_name)
    reps, lengths = orbits.representatives()
    for i, v in enumerate(reps):
        v_list = samples[i]
        assert v_list[0] == v
        g = map_leech2_vector(v)
        ref_img = triality_orbits(axis, g)
        for j, v1 in enumerate(v_list[1:]):
            if v1 == 0:
                break
            g = map_leech2_vector(v1)
            img = triality_orbits(axis, g)
            if verbose:
                print("%-3s %2d: %s" % (orbit_name, i+1, img))
            assert img == ref_img, (img, ref_img)
            if j >= N_VERIFY:
                continue
            #print("*")
            for t in T1, T2:
                assert is_good_axis(axis * g * t)
        d[(orbit_name, ref_img[0])] += lengths[i]
        d[(orbit_name, ref_img[1])] += lengths[i]
    return d


MP = True



def check_monster_axes(verbose = 0):
    orbits = load_orbits()
    samples = load_samples()
    assert isinstance(list(orbits.values())[0], Orbit_Lin2)
    print("Checking 2A axes in Monster")
    data = [(name, axis, orbits[name], samples[name])
        for name, axis in AXES.items()]
    if MP:
        with Pool() as pool:
            orders = pool.starmap(process_orbit, data)
        pool.join()
    else:
        orders = [process_orbit(*y) for y in data]
    d = defaultdict(int)
    for o in orders:
        d.update(o)

    if verbose:
        for o, n in d.items():
            print(o, n)
    with shelve.open(SHELVE_NAME) as db:
        db["mat22_suborbits"] = d
    print("Axes are as expected")    





if __name__ == "__main__":
    check_monster_axes(verbose = 0)

