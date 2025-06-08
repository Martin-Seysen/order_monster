import sys
import os
from collections import defaultdict, OrderedDict
import numpy as np
from multiprocessing import Pool
import shelve

from mmgroup import MM0, XLeech2, mat24, MM, Xsp2_Co1
from mmgroup.mm_op import mm_op_compare_abs
from mmgroup.general import Orbit_Lin2
from mmgroup.general import Orbit_Elem2


from mat22_orbits import SHELVE_NAME 
from mat22_orbits import load_orbits
sys.path.append(".")
sys.path.append(os.path.join("..", "axis_orbits"))
from utilities import compute_order


MAT22_SIZE = 22*21*20*16*3
N_TYPE42_VECTORS = 46575
CO_2_SIZE = 2**11 * MAT22_SIZE * N_TYPE42_VECTORS
H_PLUS_SIZE = 2**24 * CO_2_SIZE


def centralizer_orders(recompute = True, with_pool = True, verbose = 0):
    with shelve.open(SHELVE_NAME) as db:
        #print(list(db.keys()))
        orbit_sizes = db["ORBIT_SIZES"]

    c =  load_orbits() 
    axis_types = list(c.keys()) 
    axis_centralizers = list(c.values())
    
    del c
    d = {}
    if recompute:
        if with_pool: 
            with Pool() as pool:
                orders = pool.map(compute_order, axis_centralizers)
            pool.join()
        else:
            orders = [compute_order(x) for x in axis_centralizers]
    else:
        with shelve.open(SHELVE_NAME) as db:
            d = db["ORBIT_CENTRALIZERS"]
        orders = [d[name] for name in axis_types]

    if verbose:
        print(
"""Orders of the centralizers of the H orbits of feasible axes.

Here G_x0 is decmposed as Z = {1,x} < Q_x0 < G_x0; and we 
display the corresponding decmpositions of the centralizers."""
)    
        fmt = "%-33s %22s"
        print(fmt % ("Order of Z is:", 2))
        print(fmt % ("Order of Q_x0' = 2^{1+1+22} is:", 2**24))
        print(fmt % ("Order of Co_2 is:", CO_2_SIZE))
        print(fmt % ("Order of H = 2^{1+1+22}.Co_2 is:", H_PLUS_SIZE))
        print()
        s = "%4s %2s %9s %15s %19s   %s" % (
          "Name", "Z", "Q_x0'/Z", "H/Q_x0'", "No. of H orbits", "status"
        )
        print(s)
    for i, name in enumerate(axis_types):
        f_c, f_x, f = orders[i]
        os = orbit_sizes[name]
        product = f_c * f_x * f * os
        ok = 'ok' if product ==  H_PLUS_SIZE else 'failed'
        if verbose:
            print("%4s %2d %9d %15d %19d   %s" % 
                (name, f_c, f_x, f, os, ok)) 
        #print(product) 
        d[name] = orders[i]
    with shelve.open(SHELVE_NAME) as db:
        db["ORBIT_CENTRALIZERS"] = d
    return d


if __name__ == "__main__":
    centralizer_orders(with_pool = True)

