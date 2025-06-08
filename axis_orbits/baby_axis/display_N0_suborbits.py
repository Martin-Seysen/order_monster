"""This is a submodule for module ``watermark_suborbits``.

The main function ``check_suborbits`` checks the N_x0 orbits computed
by module ``watermark_suborbits``.
"""

import sys
import os
from collections import defaultdict, OrderedDict
import numpy as np
from argparse import ArgumentParser
import shelve

from mat22_orbits import SHELVE_NAME, load_orbits 
from watermark_suborbits import suborbit_sample_axes

sys.path.append(os.path.join("..", "utilities"))
from utilities_gap import parse_gap_output
sys.path.pop()


MAT22_SIZE = 22*21*20*16*3
N_X0_SIZE = MAT22_SIZE << 35
N_0_SIZE = (MAT22_SIZE << 34) * 6

IND_H_2B = 11707448673375
GAP_OUTPUT = os.path.join(os.path.split(__file__)[0], "Nx0_orbit_structure.txt")



def load_tables():
    global MAP_SUBORBIT, SUBORBIT_REPRESENTATIVES, SUBORBIT_SIZES
    global ORBITS, GENERATOS
    global SUBORBIT_CENTRALIZERS
    global SUBORBIT_SIZES_2
    with shelve.open(SHELVE_NAME) as db:
        MAP_SUBORBIT = db["MAP_SUBORBIT"]
        SUBORBIT_REPRESENTATIVES = db["SUBORBIT_REPRESENTATIVES"]
        SUBORBIT_SIZES = db["SUBORBIT_SIZES"]
        SUBORBIT_CENTRALIZERS = db["SUBORBIT_CENTRALIZERS"]
        SUBORBIT_SIZES_2 = db["SUBORBIT_SIZES_2"] 
        ORBITS = load_orbits()


class Nx0_Suborbit:
    """Contain elevant data of an N_x0 suborbit

    The group :math:`N_{x0}` has structure
    :math:`2^{2+11+22}.(M_{24} \times 2)`. 
    Relevant attributes of this class are:

    e:        triple of exponents corrsponding to section :math:`2^{2+11+22}`
    g_order:  order of section corresponding to :math:`M_{24}`
    s.m24_structure:
              structure of section corresponding to :math:`M_{24}`
    s_order:  exponent ot final factor group :math:`.2`
    s.images: 
              tuple of images of suborbits under the action of the
              powers of the triality element
    Gx0_orbits:
              triple of :math:`G_{x0}` corrsponding to the orbits
              given by attribute ``s.images``
    Gx0_orbits:
              triple of :math:`G_{x0}` corrsponding to the orbits
              given by attribute ``s.images``
    """
    pass



def get_Nx0_orbit(axis, e):
    wm = (axis.profile_Nxyz((e, 0), 1)[1],
           axis.profile_Nxyz((e, 1), 1)[1])
    return MAP_SUBORBIT[wm]


def key_Gx0_orbit(orbit):
    return int(orbit[:-2]) , orbit[-2], -int(orbit[-1])



def Nx0_orbit_str(axis):
    d = defaultdict(int)
    for e in range(3):
        G_x0_orbit = axis.axis_type(e)
        key = key_Gx0_orbit(G_x0_orbit)
        d[(key, G_x0_orbit, get_Nx0_orbit(axis, e))] += 1
    o_list_keyed = list(sorted(d.items()))
    o_list = [(value[1], n) for value, n in o_list_keyed]
    o_strings = [orbit + ("^%d" % n if n > 1 else "")
        for orbit, n in o_list]
    images = [value[2] for value in d] 
    o_str = ",".join(o_strings)
    o_key = (len(o_list), tuple(key_Gx0_orbit(o) for o, _ in o_list))
    return images, o_str, o_key



S_FMT_NX = "%3s %3s %3s, %3s %3s %3s, [%1s %2s %2s] %8s %1s"

def collect_Nx0_suborbits(verbose = 0):
    sub_structure = parse_gap_output(GAP_OUTPUT)
    suborbits = [None] * len(SUBORBIT_REPRESENTATIVES)
    sample_axes = suborbit_sample_axes()
    for i, _ in enumerate(SUBORBIT_REPRESENTATIVES):
        suborbits[i] = s = Nx0_Suborbit()
        gap_output = sub_structure[i]
        gap_order = gap_output.order
        s.gap_id = gap_output.id
        s.m24_structure = gap_output.structure_description()
        e, s.g_order, s.s_order = SUBORBIT_SIZES_2[i]
        s.e = [e[0] + e[1], e[2], e[3] + e[4]]
        axis = sample_axes[i]
        s.images, s.Gx0_orbits, s.Gx0_key = Nx0_orbit_str(axis)
        assert s.g_order == gap_order
        assert i in s.images
        if verbose:
            print(S_FMT_NX % tuple(s.images + s.Gx0_orbits + s.e +
               [s.g_order, s.s_order])) 
    return suborbits

class N0_orbit:
    pass





def key_N0_orbit(orbit):
    return (-orbit.order_centralizer, -orbit.g_order << sum(orbit.e),
        -sum(orbit.e), orbit.Gx0_key)

def key_m24_structure(s):
    return s.count('('), len(s), s

def find_m24_structure(suborbits, images):
    data = []
    for i in set(images):
        s = suborbits[i]
        data.append((i, s.m24_structure, s.gap_id, s.g_order))
    ind, struct, id, ord = list(zip(*data))
    neq = sum(struct[i-1] != struct[1] for i in range(1, len(struct)))
    if not neq:
        return struct[0]
    struct = sorted(struct, key = key_m24_structure)
    neq_orders = sum(ord[i-1] != ord[1] for i in range(1, len(ord)))
    if neq_orders:
        o = dict(zip(ind, ord))
        ERR = "Groups have different orders: %s"
        raise ValueError(ERR % o)
    neq = sum(id[i-1] != id[1] for i in range(1, len(id)))
    if neq or id[0] == None:
        W = "Warning: group structures might be unequal (order = %s)"
        print(W % ord[0]) 
        print(" ", struct)
    return struct[0]
        
    


def join_Nx0_orbits(suborbits, verbose = 0):
    orbits = []
    used = set()
    num_axes = 0
    num_Nx0_orbits = 0
    for i, suborbit in enumerate(suborbits):
        if i in used:
            continue
        used.add(i)
        orb = N0_orbit()
        orbits.append(orb)
        n = len(orbits)
        orb.Gx0_orbits = suborbit.Gx0_orbits
        orb.Gx0_key = suborbit.Gx0_key
        orb.e = suborbit.e
        orb.g_order = suborbit.g_order
        orb.m24_structure = find_m24_structure(suborbits, suborbit.images)
        num_S3_cosets = 0
        orb.Nx0_orbits = suborbit.images
        for j in set(suborbit.images):
            used.add(j)
            img_suborbit = suborbits[j]
            assert orb.Gx0_orbits == suborbit.Gx0_orbits
            assert orb.e == img_suborbit.e
            assert orb.g_order == img_suborbit.g_order
            num_S3_cosets += 2 >> img_suborbit.s_order
        num_Nx0_orbits += len(suborbit.images)
        assert  num_S3_cosets in [1,2,3,6]
        orb.s_order = 6 // num_S3_cosets
        orb.order_centralizer = (orb.g_order * orb.s_order) << sum(orb.e)
        quot, rem = divmod(N_0_SIZE, orb.order_centralizer)
        assert rem == 0
        num_axes += quot
        if verbose > 1:
            print(n, orb.Gx0_orbits, orb.e, orb.g_order, orb.s_order,
                orb.m24_structure)
    assert num_axes == IND_H_2B, num_axes
    assert num_Nx0_orbits == len(suborbits)
    orbits = sorted(orbits, key = key_N0_orbit)
    if verbose:
        for i, orb in enumerate(orbits):
            n = i + 1
            print(n, orb.Gx0_orbits, orb.e, orb.g_order, orb.s_order,
               orb.m24_structure)
    return orbits        
           



def format_Nx0_orbits(suborbits, N0_orbits):
    print("""
Decomposition of N_x0 orbits of 2A axes into N_xyz orbits

Orbit-Number  Size of C_x0       Index C_x0/C_xyz    G_x0
N_x0  N_0      e        |G|      ind  2^e   G        orbit""")
    FMT = "%4d %4d    %3d %10d       %1d    %1d    %1d        %-5s"
    def index2(a, b):
        index, rem = divmod(a, b)
        assert rem == 0, (a, b, index, rem)
        assert 1 <= index <= 2, (a, b, index, rem)
        return index
    from mmgroup import Xsp2_Co1
    from utilities import order_Nx0
    sample_axes = suborbit_sample_axes()
    map_N0 = {}
    for j, orb in enumerate(N0_orbits):
        for  Nx0_orbit in orb.Nx0_orbits:
            map_N0[Nx0_orbit] = j + 1, orb
    for i, _ in enumerate(SUBORBIT_REPRESENTATIVES):
        e, g, s = SUBORBIT_SIZES_2[i]
        index = 1 << s
        j, orb = map_N0[i]
        assert sum(e) == sum(orb.e)
        assert g == orb.g_order
        Nxyz_e = sum(orb.e)
        Nxyz_g_order = orb.g_order
        c = [Xsp2_Co1(g) for g in SUBORBIT_CENTRALIZERS[i]]
        e, Nx0_g_order = order_Nx0(c , e_odd = False)
        Nx0_e = sum(e)
        index_Nxyz_e = 1 << (Nx0_e - Nxyz_e)
        assert index_Nxyz_e in [1, 2]
        index_Nxyz_g = index2(Nx0_g_order, Nxyz_g_order) 
        assert index_Nxyz_g * index_Nxyz_e == index
        axis = sample_axes[i]
        Gx0_orbit = axis.axis_type()
        print(FMT % (i, j, Nx0_e, Nx0_g_order,
             index, index_Nxyz_e, index_Nxyz_g, Gx0_orbit ))


def format_N0_orbits(orbits, tex, f=sys.stdout):
    HD = """Table of N_0 orbits of 2A axes

 No  e            G                                     |G|  S    G_x0"""
    FMT = "%3d  %-12s %-32s %8s  %-3s  %-4s" 
    if not tex:
        print(HD) 
    else:
        print("Table of N_0 orbits of 2A axes\n")
    imax = len(orbits) - 1
    for i, orb in enumerate(orbits):
        n = i+1
        e = orb.e[:]
        while len(e) and e[0] == 0:
            e = e[1:]
        if len(e):
            if len(e) == 1 and e[0] == 1:
                e_str = "2"
            else:
                e_str = "2^{%s}" % "+".join(map(str, e)) 
        else:
            e_str = "1"
        g_str = orb.m24_structure
        ord_g_str = str(orb.g_order)
        s_str = "S_3" if orb.s_order == 6 else str(orb.s_order)
        Gx0_str = orb.Gx0_orbits
        if tex:
            s =  f"${n}$ & ${e_str}$ & ${g_str}$ & ${ord_g_str}$ &"
            s += f" ${s_str}$ & ${Gx0_str}$"
            if  i != imax: 
                s += r"\\"
        else:
            s = FMT % (n, e_str, g_str, ord_g_str, s_str, Gx0_str)
        print(s, file = f)



def parse_args():
    description = ('Display N_0 orbits of 2A axes. ')
    usage = "usage: %prog [options]"
    parser = ArgumentParser(usage, description = description)
    parser.add_argument("-t",  dest="latex", action="store_true",
        help = "Output data in format suitable for LaTex")
    parser.add_argument("-x",  dest="N_x0", action="store_true",
        help = "Display mapping from N_x0 orbits to N_0 orbits")
    options  = parser.parse_args()
    return options



SUBORBITS = None
N0_ORBITS = None


def load_suborbits():
    global SUBORBITS, N0_ORBITS
    load_tables()
    if SUBORBITS is None:
        SUBORBITS = collect_Nx0_suborbits()
    if N0_ORBITS is None:
        N0_ORBITS = join_Nx0_orbits(SUBORBITS)
    return SUBORBITS, N0_ORBITS


def display_N0_orbits(latex):
    _, N0_ORBITS = load_suborbits()
    format_N0_orbits(N0_ORBITS, latex)

def display_Nxyz_orbits(latex):
    SUBORBITS, N0_ORBITS = load_suborbits()
    format_Nx0_orbits(SUBORBITS, N0_ORBITS)



if __name__ == "__main__":
    opt = parse_args()
    display_N0_orbits(opt.latex)
    if opt.N_x0:
        print(
        "----------------------------------------------------------------"
        )
        display_Nxyz_orbits(opt.latex)


