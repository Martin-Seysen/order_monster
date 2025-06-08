import sys
import shelve
from collections import defaultdict, OrderedDict
import numpy as np
from argparse import ArgumentParser

import sympy
from sympy import *

sys.path.append(r".")

from mmgroup import MM0
from mmgroup.axes import Axis
from mat24_orbits import SHELVE_NAME 
from mat24_orbits import AXES 
from mat24_orbits import load_orbits
ORBITS = list(AXES.keys())




DICT_NAME = "mat24_suborbits"
AXES = Axis.representatives()



ORBITS_NORTON = {
 '2A': 2**4 * 3**3 * 5 * 7 * 13,
 '2B': 2**8 * 3**4 * 5**2 * 7 * 11 * 13 * 23,
 '4A': 2**16 * 3**7 * 5**3 * 7 * 13,
 '4B': 2**21 * 3**5 * 5**3 * 7 * 11 * 13 * 23,
 '4C': 2**20 * 3**7 * 5**3 * 7 * 11 * 13 * 23,
 '6A': 2**28 * 3**3 * 5**3 * 7 * 13 * 23,
 '6C': 2**28 * 3**4 * 5**3 * 7**2 * 11 * 13 * 23,
 '6F': 2**31 * 3**5 * 5**3 * 7 * 11 * 13 * 23,
 '8B': 2**31 * 3**7 * 5**3 * 7**2 *  13 * 23,
'10A': 2**35 * 3**7 * 5 * 7 *  13 * 23,
'10B': 2**32 * 3**7 * 5**2 * 7**2 * 11 * 13 * 23,
'12C': 2**36 * 3**5 * 5**3 * 7 * 11 * 13 * 23,
}


GROUPS_CO_1_MAGMA_TEX = {
 '2A':  r'\mbox{Co}_2',
 '2B':  r'2^{1+8}.O_8^+(2)',
 '4A': r'2^{11}.M_{23}',
 '4B': r'2^{1+8}.S_6(2)',
 '4C': r'2^{1+8+6}.A_8',
 '6A': r'U_6(2).2',
 '6C': r'2^{1+8}.(3 \times U_4(2)).2',
 '6F': r'2^{1+8}.A_9',
 '8B': '2^{1+10}.M_{11}',
'10A':  r'\mbox{HS}.2',
'10B': r'2^{1+8}.(A_5 \times A_5).2',
'12C': r'2.S_6(2)',
}

with shelve.open(SHELVE_NAME) as db:
    monster_tables = db[DICT_NAME]
    orbits = ORBITS
    #print(orbits)
    #print(monster_tables)
    m = np.zeros((12,12), dtype = np.uint32)
    for i in range(12):
        s = 0
        for j in range(12):
            m[i,j] = x =  monster_tables[(orbits[i],orbits[j])]
            s  += m[i,j]
        assert s  == 16584750, (i, j, s)
    m1 = m - 16584750 * eye(12) 


def involution_type(axis):
    inv = axis.central_involution()
    inv_type, _ = inv.conjugate_involution_G_x0()
    group = axis.auto_group
    powers = axis.powers
    pos = powers.rfind(',')
    powers = powers[:pos+1] + inv_type
    return (powers, group)


IND_2B_M = 97239461142009186000




def compute_orbits(verbose = 0):
    #print(m)
    ns = m1.T.nullspace()
    assert len(ns) == 1
    f = 196560 / ns[0][0]
    #print("FFFF", f, type(ns[0]))
    nsf = ns[0]*f
    s = 0
    orbit_sizes = {}
    if verbose:
        print("%4s: %21s   %-12s  %-22s" % (
            'type', 'orbit size', 'powers', 'centralizer'))
    for i in range(12):
        x = nsf[i]
        s += x
        axis = AXES[ORBITS[i]]
        powers, group = involution_type(axis)
        if verbose:
            print("%4s: %21d   %-12s  %-22s" % (
                ORBITS[i], x, powers, group))
        if ORBITS[i] in ORBITS_NORTON:
            xn =  ORBITS_NORTON[ ORBITS[i] ]
            assert x == xn, xn
        else:
            print(" ?")
        orbit_sizes[ORBITS[i]] = x
    if verbose:
        print("\n sum: %21d" % s)
    assert (m.T @ nsf == 16584750 * nsf)
    if verbose > 1:
        print("The relevant eigenvector of the matrix has shape",
              nsf.shape)
    with shelve.open(SHELVE_NAME) as db:
        db["ORBIT_SIZES"] = orbit_sizes
    assert s == IND_2B_M
    return s



def load_centralizers():
    try:
        with shelve.open(SHELVE_NAME) as db:
            return db["ORBIT_CENTRALIZERS"]
    except:
        from centralizer_orders import centralizer_orders
        centralizer_orders(with_pool = False, verbose = 0)
        with shelve.open(SHELVE_NAME) as db:
            return db["ORBIT_CENTRALIZERS"]

def load_orbit_sizes():
    with shelve.open(SHELVE_NAME) as db:
        return db["ORBIT_SIZES"]

HD = "Sizes and centralizers of the orbits of G_x0 on the 2A axes\n"


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

 
def display_orbits(Nxyz = False):
    print(HD)
    s = s_o = s_xyz = 0
    orbit_sizes = load_orbit_sizes()
    centralizers = load_centralizers()
    orbits = load_orbits()
    hd = ('type', '#N_x0', 'orbit size', 'centralizer', 'powers')
    if Nxyz:
        FMT = "%4s: %4s %4s  %20s  %-22s  %-12s"
        FMT_HD = "%4s: %9s  %19s  %-22s  %-12s"
        hd1 = (hd[0], 'No. orbits') + hd[2:]
        d_xyz = load_num_N_xyz_suborbits()
        hd = hd[:2] + ('#N_xyz',) + hd[2:]
        print(FMT_HD % hd1)
        print( "%4s  %4s %4s" % ("", "N_x0", "N_xyz"))
    else:
        FMT = "%4s: %5s  %20s  %-22s  %-12s"
        print(FMT % hd)
    for i, name in enumerate(ORBITS):
        x = orbit_sizes[name]
        s += x
        n_orbits = orbits[name].n_orbits()
        s_o += n_orbits
        axis = AXES[ORBITS[i]]
        powers, group = involution_type(axis)
        data =  (name, n_orbits, x, group, powers)
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
                GROUPS_CO_1_MAGMA_TEX[name], 
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
   compute_orbits(verbose = 0)
   if options.tex:
       display_orbits_tex()
   elif options.pydict:
       display_orbits_pydict()
   else:
       display_orbits()









