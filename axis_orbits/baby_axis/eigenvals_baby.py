import sys
import shelve
import numpy as np
from collections import defaultdict, OrderedDict
from argparse import ArgumentParser

import sympy
from sympy import *

sys.path.append(r".")

from mmgroup import MM0
from mmgroup.axes import BabyAxis
from mat22_orbits import SHELVE_NAME 
from mat22_orbits import AXES 
from mat22_orbits import load_orbits
ORBITS = list(AXES.keys())



DICT_NAME = "mat22_suborbits"
AXES = BabyAxis.representatives()


GROUPS_H_MAGMA_TEX = {
 '2A1':  r'\mbox{Co}_2',
 '2A0':  r'2^{10}.M_{22}.2',
 '2B1':  r'2^{1+8}.S_6(2)',
 '2B0':  r'2^{1+4+6}.A_8',
 '4A1': r'2^{9+1}.L_3(4).2',
 '4B1': r'2^{1+8+5}.S_6',
 '4C1': r'2^{1+4+6}.A_8',
 '6A1': r'U_6(2).2',
 '6C1': r'2^{1+8}.U_4(2).2',
'10A1':  r'\mbox{HS}.2',
}


with shelve.open(SHELVE_NAME) as db:
    monster_tables = db[DICT_NAME]
    orbits = ORBITS
    #print(orbits)
    #print(monster_tables)
    m = np.zeros((10,10), dtype = np.uint32)
    for i in range(10):
        s = 0
        for j in range(10):
            m[i,j] = x =  monster_tables[(orbits[i],orbits[j])]
            s  += m[i,j]
        assert s  == 93150, (i,s)
    m1 = m - s * eye(10) 


def compute_orbits(verbose = 0):
    S = "Sizes of orbits of 2A axes orthogonal to v^+ under the action of H"
    if verbose:
        print("\n%s:\n" % S)
    #print(m)
    ns = m1.T.nullspace()
    assert len(ns) == 1
    f =  1 / ns[0][0]
    #print("FFFF", f, type(ns[0]))
    nsf = ns[0]*f
    s = 0
    orbit_sizes = {}
    if verbose:
        print("%4s: %15s" % ('type', 'orbit size'))
    for i in range(10):
        x = nsf[i]
        s += x
        axis = AXES[ORBITS[i]]
        #powers, group = involution_type(axis)
        if verbose:
            print("%4s: %15d" % (ORBITS[i], x))
        orbit_sizes[ORBITS[i]] = x
    if verbose:
        print("")
        print(" sum: %15d" % s)

    with shelve.open(SHELVE_NAME) as db:
        db["ORBIT_SIZES"] = orbit_sizes

    IND_G_2B = 11707448673375

    assert s == IND_G_2B




def load_centralizers():
    try:
        with shelve.open(SHELVE_NAME) as db:
            return db["ORBIT_CENTRALIZERS"]
    except:
        from centralizer_orders import centralizer_orders
        centralizer_orders(with_pool = False)
        with shelve.open(SHELVE_NAME) as db:
            return db["ORBIT_CENTRALIZERS"]

def load_orbit_sizes():
    with shelve.open(SHELVE_NAME) as db:
        return db["ORBIT_SIZES"]


HD = "Sizes and centralizers of the orbits of H on feasible 2A axes\n"


def load_num_N_xyz_suborbits():
    """Count N_xyz suborbits in a G_x0 orbit

    The function returns a dictionary mappeng the G_x0 orbits
    to the number of N_xyz_suborbit contained in that orbit.

    Caution:

    Calling this function requires the shelve entries SUBORBIT_SIZES_2
    and SUBORBIT_REPRESENTATIVES. These entries are computed after
    calling the main function ``compute_orbits`` of thie module!
    """
    d = defaultdict(int)
    with shelve.open(SHELVE_NAME) as db:
        representatives = db["SUBORBIT_REPRESENTATIVES"]
        sizes_2 = db["SUBORBIT_SIZES_2"]
    for i, (orbit_name, entry, v) in enumerate(representatives):
        axis = AXES[orbit_name] * MM0('c', v) ** -1
        G_x0_orbit = axis.axis_type()
        n = 2 >> sizes_2[i][2]
        d[G_x0_orbit] += n
    return d





def str_Q_x0(orders):
    l0 = orders[0].bit_length()-1
    l1 = orders[1].bit_length()-1
    if l0  == 0:
        if l1  <= 1:
            return "$%d$" % (2**l1)
        else:
            return "$2^{%d}$" % l1
    else:
       return "$2^{%d+%d}$" % (l0, l1)


def display_orbits(Nxyz = False):
    print(HD)
    s = s_o = s_xyz = 0
    orbit_sizes = load_orbit_sizes()
    centralizers = load_centralizers()
    orbits = load_orbits()
    hd = ('type', '#N_x0', 'orbit size', 'centralizer')
    if Nxyz:
        FMT = "%4s: %4s %4s  %20s  %-22s"
        FMT_HD = "%4s: %9s  %19s  %-22s"
        hd1 = (hd[0], 'No. orbits') + hd[2:]
        d_xyz = load_num_N_xyz_suborbits()
        hd = hd[:2] + ('#N_xyz',) + hd[2:]
        print(FMT_HD % hd1)
        print( "%4s  %4s %4s" % ("", "N_x0", "N_xyz"))
    else:
        FMT = "%4s: %5s  %20s  %-22s"
        print(FMT % hd)
    for i, name in enumerate(ORBITS):
        x = orbit_sizes[name]
        s += x
        n_orbits = orbits[name].n_orbits()
        s_o += n_orbits
        axis = AXES[ORBITS[i]]
        group = str_Q_x0(centralizers[name])[1:-1]
        group += "." +  GROUPS_H_MAGMA_TEX[name]
        data =  (name, n_orbits, x, group)
        if Nxyz:
            n_xyz = d_xyz[name]
            s_xyz += n_xyz
            data = data[:2] + (n_xyz,) +  data[2:] 
        print(FMT % data)
    if Nxyz:
        print("\n sum: %4s %4s  %20s" % (s_o, s_xyz, s))
    else:
        print("\n sum: %5s  %20s" % (s_o, s))



def display_orbits_tex(Nxyz = False):
    print(HD)
    s = s_o = s_xyz = 0
    orbit_sizes = load_orbit_sizes()
    centralizers = load_centralizers()
    orbits = load_orbits()
    if Nxyz:
        d_xyz = load_num_N_xyz_suborbits()
    print("Table content for LaTeX")
    total_orbits = 0
    for i, name in enumerate(ORBITS):
        x = orbit_sizes[name]
        s += x
        n_orbits = orbits[name].n_orbits()
        s_o += n_orbits
        data = (name, n_orbits, orbit_sizes[name], 
                str_Q_x0(centralizers[name]),
                GROUPS_H_MAGMA_TEX[name], 
            )
        if Nxyz:
            n_xyz = d_xyz[name]
            s_xyz += n_xyz
            data = data[:2] + (n_xyz,) +  data[2:] 
            print("%s & %d & %d & %d &  %s  & $%s$  \\\\" % data)
        else:
            print("%s & %d & %d &  %s  & $%s$  \\\\" % data)
    print(r"\hline")
    if Nxyz:
       print(r"Total & %d & %d & %d \\" % (s_o,s_xyz, s))   
    else: 
       print(r"Total & %d & %d\\" % (s_o, s))   




def display_orbits_pydict():
    with shelve.open(SHELVE_NAME) as db:
        orbit_sizes = db["ORBIT_SIZES"]
    print("{")
    for name, size in orbit_sizes.items():
        print("'%s' : %d," % (name,size))
    print("}")


def parse_args():
    description = ("""Display G_x0 orbits on 2A axes in the Monster"""
    )
    usage = "usage: %prog [options]"
    parser = ArgumentParser(usage, description = description)
    parser.add_argument("-t",  dest="tex", action="store_true",
        help = "Output data also in LaTex format")
    parser.add_argument("-d",  dest="pydict", action="store_true",
        help="Output a python dintionary containing orbit sizes")
    parser.add_argument("-r",  dest="run", action="store_true",
        help = "Explain prerequisites for running this script and exit")    
    options  = parser.parse_args()
    return options



def explain_prerequisites():
    print("""You must run the following python scripts before this script:
python mat24_orbits.py
python suborbits.py
""")



def explain_prerequisites():
    print("""You must run the following python scripts before this script:
python mat22_orbits.py
python suborbits.py
""")



def show_eigenvals(latex, Nxyz = False):
   if latex:
       display_orbits_tex(Nxyz)
   else:
       display_orbits(Nxyz)



if __name__ == "__main__":
   options = parse_args()
   if options.run:
       explain_prerequisites()
       sys.exit(0)
   compute_orbits()
   if options.tex:
       display_orbits_tex()
   if options.pydict:
       display_orbits_pydict()











