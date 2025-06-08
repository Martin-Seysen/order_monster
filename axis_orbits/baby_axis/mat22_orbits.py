import os
import time
import shutil
import numpy as np
from argparse import ArgumentParser
import shelve
from multiprocessing import Pool

from mmgroup import MM0, XLeech2, leech2_orbits_raw, mat24
from mmgroup.axes import  BabyAxis, set_axis_group

from mmgroup import MM0, XLeech2, leech2_orbits_raw, mat24, Xsp2_Co1
from mmgroup.axes import  Axis, set_axis_group
from mmgroup.general import Orbit_Lin2

try:
    import mmgroup.mm_reduce
    from mmgroup import MM
except (ImportError, ModuleNotFoundError):
    MM = MM0

def configure_axis_group():
    set_axis_group(group = MM0, shorten = False) 

configure_axis_group()


#set_axis_group(group = MM0, shorten = True) 
AXES = BabyAxis.representatives()
G = BabyAxis.group
BETA_PLUS = XLeech2(0x200)

orbits = list(AXES.keys())
#print(orbits)


SHELVE_PATH = os.path.join(os.path.split(__file__)[0], "shelve")
SHELVE_NAME = os.path.join(SHELVE_PATH,"mat22_orbit_tables")

##################################################################




def find_axis_centralizer(orbit, verbose = False, g_start = None):
    g1 = g_start if g_start is not None else  G('r', 'G_x0 & B')
    rep = AXES[orbit].copy()
    ax = rep * g1
    g2 = ax.reduce_G_x0()
    assert rep * g1 * g2 == rep
    return g1 * g2

def display_general_header():
    s = r"""
Generators of centralizers of 2A axes in H = G_x0 \cap H^+.
Here H^+ is the centralizer of the 2A axis v^+.
We consider 2A axes orthogonal to the axis v^+ only.
We display one reprpresentative of each G_x0 orbit of 2A axes.
Central involution z of H: %s
Standard 2A involution \beta^- for the 2A axis v^- : %s
"""
    std_ax = AXES['2A1']
    print(s % (std_ax.g_central, std_ax.g_axis_start))



def display_orbit_header(orbit):
    s = """
Orbit %s
Representative of orbit is t = \\beta^- ** h
t = %s
h = %s
Central involution i = (t * z) ** %d of dihedral group <t, z>
i = %s
mmgroup character chi_G_x0 of involution i: %s
class of involution i: %s
"""
    ax = AXES[orbit]
    z = MM(ax.g_central)
    h = MM(ax.g)
    t = MM(ax.g_axis)
    o, i = (t * z).half_order()
    ch = i.chi_G_x0()
    iclass, _ = i.conjugate_involution_G_x0()
    print(s % (orbit, t, h, o//2, i, ch, iclass))



def make_generators(n_generators = 20, verbose = 0):
    #initialize_all()
    if verbose:
        print("Orbits analysed:\n%s" % orbits)
    d = {}
    for orbit in orbits:
        if verbose:
            for i in range(3):
                 AXES[orbit].display_sym(i, text="t**%d" % i)
            print(orbit)
        gen = []
        for i in range(n_generators):
            gen.append(find_axis_centralizer(orbit, verbose))
        if not None in gen:
            d[orbit] = gen
    if verbose:
        print("Generators of centralizers found for orbits:")
        print(list(d.keys()))
    return d

########################################################################


def _pickle_generators(generator_list):
    a = np.zeros((len(generator_list), 10), dtype = np.uint32)
    for i, gen in enumerate(generator_list):
        mm = Xsp2_Co1(gen).mmdata
        a[i, :len(mm)] = mm
    return a

def _unpickle_generators(a):
    return [Xsp2_Co1('a', data) for data in a]

_PIC = _pickle_generators, _unpickle_generators

def _map_generator(g):
    return g.as_compressed_Co1_bitmatrix()


MAX_ORBIT_SAMPLE_SIZE = 8

def type42_orbits_samples(orbits):
    BETA = BETA_PLUS.ord
    assert isinstance(orbits, Orbit_Lin2)
    reps = [v for v in orbits.representatives()[0]
        if XLeech2(v).type == 4 and XLeech2(v ^ BETA).type == 2]
    ls = MAX_ORBIT_SAMPLE_SIZE + 2
    orbit_samples = np.zeros((len(reps), ls), dtype = np.uint32)
    for i, v in enumerate(reps):
        orbit = orbits.orbit(v)
        size = min(len(orbit), MAX_ORBIT_SAMPLE_SIZE)
        sample = np.random.choice(orbit, size, replace = False) 
        orbit_samples[i, 0] = v
        orbit_samples[i, 1 : len(sample) + 1] = sample
    return orbit_samples


def get_orbits(generators):
    orbits = Orbit_Lin2(_map_generator, generators)
    orbit_samples = type42_orbits_samples(orbits)
    compressed_orbits = orbits.compress(orbit_samples[:, 0])
    return compressed_orbits, orbit_samples

PICKLE_FUNTIONS = None

def store_pickle_functions():
    global PICKLE_FUNTIONS
    if PICKLE_FUNTIONS is None: 
        orbits = Orbit_Lin2(_map_generator)
        _, PICKLE_FUNTIONS = orbits.pickle(*_PIC)

MP = True


def compute_orbits(n_generators = 10, store = True, verbose = 0):
    store_pickle_functions() 
    d = make_generators(n_generators, verbose)
    list_orbit_names, list_generators = list(d.keys()), list(d.values()) 
    if MP:
        with Pool() as pool:
            orbits_samples = pool.map(get_orbits, list(list_generators))
        pool.join()
    else:
        orbits_samples = [get_orbits(y) for y in list_generators]
    list_orbits = [orbits for orbits, _ in orbits_samples]
    list_samples = [samples for _, samples in orbits_samples]
    d = dict(zip(list_orbit_names, list_orbits))
    d_samples = dict(zip(list_orbit_names, list_samples))
    if store:
         d_pic = {}
         for name, obj in d.items():
             data, _ = obj.pickle(*_PIC)
             d_pic[name] = data
         with shelve.open(SHELVE_NAME) as db:
             db["Lin2Orbits"] = d_pic
             db["Lin2Samples"] = d_samples
         del d_pic
    return d    

def load_orbits():
    store_pickle_functions()
    with shelve.open(SHELVE_NAME) as db:
        d0 = db["Lin2Orbits"]
    d = {}
    for name, pickled in d0.items():
        d[name] = Orbit_Lin2(pickled, PICKLE_FUNTIONS)
    del d0
    return d


def load_samples():
    store_pickle_functions()
    with shelve.open(SHELVE_NAME) as db:
        d = db["Lin2Samples"]
    return d



########################################################################




def parse_args():
    description = ('Generate centralizers of representatives of H '
    'orbits of feasible 2A axes. ' 
    )
    usage = "usage: %prog [options]"
    parser = ArgumentParser(usage, description = description)
    parser.add_argument("-n",  dest="no_mat24", action="store_true",
        help = "Do not display information about orbits in Mat24")
    parser.add_argument("-r",  dest="recompute", action="store_true",
        help="Recompute all precomputed data")
    parser.add_argument("-s",  dest="sizes", action="store_true",
        help = "Display sizes of orbits in Mat24")
    parser.add_argument("-v",  dest="verbose", action="store_true",
        help="Verbose operation" )
    
    options  = parser.parse_args()
    return options



def check_recompute(recompute = False):
    if recompute and os.path.exists(SHELVE_PATH):
        shutil.rmtree(SHELVE_PATH)   # Remove directory if it exists
    if not os.path.exists(SHELVE_PATH):
        os.makedirs(SHELVE_PATH)
        recompute = True
    return recompute


def display_orbits(d, header = True, mat24 = True, sizes = False):
    display_general_header()
    for orbit_name, axis in AXES.items():
        display_orbit_header(orbit_name)
        orbits = d[orbit_name]  
        if mat24:
            S = "G_x0 orbit %-3s contains %4d N_x0 orbits"
            print(S % (orbit_name, orbits.n_orbits()))
            if sizes:
                orbit_sizes = list(orbits.representatives()[1])
                orbit_sizes.sort()         
                orbit_sizes.reverse()         
                print("Orbit sizes")
                for i, s in enumerate(orbit_sizes):
                    print("%2d:%8d" % (i+1, s))
        print("Generators for orbit centralizer:")
        v = axis.v15
        for g in orbits.generators():
            print(MM(g))

def check_orbits(d, verbose = 0):
    for orbit_name, axis in AXES.items():
        orbits = d[orbit_name] 
        v = axis.v15
        for g in orbits.generators():
            assert v * G(g) == v, (orbit, g,  (v * G(g))['B',:4,:4], v['B',:4,:4])
    if verbose:
        print("Centralizers of axes are correct")


if __name__ == "__main__":
    options = parse_args()
    recompute = check_recompute(options.recompute)
    if recompute:
        d = compute_orbits(n_generators = 10, store = True) 
    d = load_orbits()   
    check_orbits(d)
    display_orbits(d, not options.no_mat24, options.sizes)


