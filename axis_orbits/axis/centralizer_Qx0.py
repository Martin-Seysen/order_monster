import sys
import os
from collections import defaultdict, OrderedDict
import numpy as np

from mmgroup import MM0, XLeech2, leech2_orbits_raw, mat24, MM
from mmgroup.mm_op import mm_op_compare_abs
from mmgroup.axes import Axis



AXES = Axis.representatives()





def centralize_axis_Qx0(axis):
    v15 = axis.v15
    ax0 = [v15 * MM(XLeech2(i)) for i in range(1 << 12)]
    ax1 = [v15 * MM(XLeech2(i)) for i in range(0, 1 << 25, 1 << 12)]
    hashes = defaultdict(lambda : ([], []))
    for ax in ax0:
        h = hashes[ax.hash()]
        h[0].append(ax)
    #print(".", end="", flush=True)
    for ax in ax1:
        h = hashes[ax.hash()]
        h[1].append(ax)
    #print(".", end="", flush=True)
    for h, (a0, a1) in hashes.items():
        if len(a0)  * len(a1):
            # print(len(a0), len(a1))
            pass
        if len(a0) * len(a1) > 1 << 18:
            return None
    #print(".", end="", flush=True)
    n = 0
    for h, (a0, a1) in hashes.items():
        #print(len(a0), len(a1))
        for v0 in a0:
            for v1 in a1: 
                if v0 == v1:
                    n += 1
    #print(".", end="", flush=True)
    return n
        


for orbit, ax in AXES.items():
    #print(orbit, ax.__dict__.keys())
    n = centralize_axis_Qx0(ax)
    # print(orbit)
    if n is not None:
        print("%-3s: %6d" % (orbit, n))

   
 
        
    






