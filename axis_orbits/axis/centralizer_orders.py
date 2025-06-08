import sys
import os
from collections import defaultdict, OrderedDict
import numpy as np
from multiprocessing import Pool
import shelve

from mmgroup import MM0, XLeech2, mat24, MM, Xsp2_Co1


from mat24_orbits import SHELVE_NAME 
from mat24_orbits import load_orbits

sys.path.append(os.path.join("..", "utilities"))
from utilities import order_Nx0, compute_order
sys.path.pop()


MAT24_SIZE = 24*23*22*21*20*16*3
N_TYPE4_VECTORS = 8292375
CO_1_SIZE = 2**11 * MAT24_SIZE * N_TYPE4_VECTORS
G_X0_SIZE = 2**25 * CO_1_SIZE


        

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
"""Orders of the centralizers of the G_x0 orbits of the axes.

Here G_x0 is decmposed as Z = {1,x} < Q_x0 < G_x0; and we 
display the corresponding decmpositions of the centralizers."""
) 
        print("Order of G_x0 is:", G_X0_SIZE)
        print()
        s = "%4s %2s %8s %15s %21s   %s" % ("Name", "Z", "Q_x0/Z",
        "G_x0/Q_x0", "No. of G_x0 orbits", "status"
        )
        print(s)
        
    for i, name in enumerate(axis_types):
        f_c, f_x, f = orders[i]
        os = orbit_sizes[name]
        product = f_c * f_x * f * os
        ok = 'ok' if product ==  G_X0_SIZE else 'failed'
        if verbose:
            print("%4s %2d %8d %15d %21d   %s" % 
                (name, f_c, f_x, f, os, ok)) 
        #print(product) 
        d[name] = orders[i]
    if recompute:
        with shelve.open(SHELVE_NAME) as db:
            db["ORBIT_CENTRALIZERS"] = d
    return d


if __name__ == "__main__":
    centralizer_orders(recompute = True, verbose=1)

