import sys
import os
import time
from collections import defaultdict, OrderedDict
from random import randint, sample
import numpy as np
from multiprocessing import Pool
import shelve



from mmgroup import MM0, XLeech2, leech2_orbits_raw, mat24, MM
from mmgroup.general import Orbit_Lin2

from mat24_orbits import AXES, Axis, configure_axis_group
from mat24_orbits import SHELVE_NAME 
from mat24_orbits import load_orbits, load_samples
from mat24_orbits import _map_generator




def make_axis_orbit(axis, orbits):
    assert isinstance(axis, Axis)
    assert isinstance(orbits, Orbit_Lin2)
    generators = orbits.generators()
    sizes = sorted(list(orbits.representatives()[1]))
    #print(sizes)
    n_generators = 2
    for i in range(20):
         new_generators = sample(generators, 
                             min(n_generators, len(generators)))
         new_orbits = Orbit_Lin2(_map_generator, new_generators)
         new_reps, new_sizes = new_orbits.representatives()
         new_data = [(v, n) for v, n in zip(new_reps, new_sizes)
               if XLeech2(v).type == 4]
         new_sizes = sorted([x[1] for x in new_data])
         if sorted(list(new_sizes)) ==  sizes:
              new_reps = [x[0] for x in new_data]
              new_orbits.compress(new_reps)
              return new_orbits, new_reps
         if i == 4:
             n_generators += 1
  
    err = "No small set of generators for centralizer of axis found"
    raise ValueError(err)   

     



def centralizers(orbit_name):
    return orbits[orbit_name].generators()



def map_leech2_vector(v):
    return MM0('c', v) ** -1


def mmstr(g):
    return str(MM(g))[1:]


T1 = MM0('t', 1); T2 = MM0('t', 2) 


def make_certificate_axis(name, axis, axis_orbit):
    cert = []
    cent_strings = set()
    orbits, reps = make_axis_orbit(axis, axis_orbit)
    gs = mmstr(axis.g)
    cert.append(f"axis: {name} {gs}")
    for g in orbits.generators():
         gs = mmstr(g)
         cert.append(f"cent: 1 {gs}")
         cent_strings.add(gs)
    for g in axis_orbit.generators():
         gs = mmstr(g)
         if gs not in cent_strings:      
             cert.append(f"cent: 2 {gs}")
    for v in reps:
         g = map_leech2_vector(v) 
         gs = mmstr(g)
         orbit_size = orbits.orbit_size(v)
         cert.append(f"orb:  {orbit_size}  {gs}")
         ax = axis * g
         for exp, t in zip([1,2], [T1, T2]):
             ax_t = ax * t
             ax_t_type = ax_t.axis_type()
             h = ax_t.reduce_G_x0()
             hs = mmstr(h)
             cert.append(f"tau{exp}: {ax_t_type}  {hs}")
    return "\n".join(cert) + "\nend:\n"



def compute_certificate():
    configure_axis_group() 
    all_orbits = load_orbits()
    pool_data = []
    for name, axis in AXES.items():
        orbit = all_orbits[name]
        pool_data.append((name, axis, orbit))

    with Pool() as pool:
       cert_list = pool.starmap(make_certificate_axis, pool_data)
    pool.join()
    cert = "".join(cert_list)
    return cert


def make_certificate(certificate_path, verbose = 0):
    s = "Computing a certificate for checking the order of the Monster"
    print(s + "...")
    cert = compute_certificate()
    if verbose:
        print(cert)
    with open(certificate_path, "wt") as f:
        f.write(cert)
    print("Certificate computed")


