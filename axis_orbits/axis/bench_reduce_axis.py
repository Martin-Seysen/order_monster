import sys
import os
from collections import defaultdict, OrderedDict
import numpy as np
import time

from mmgroup import MM0, XLeech2, leech2_orbits_raw, mat24, MM
from mmgroup.mm_op import mm_op_compare_abs
from mmgroup.axes import Axis



STD_AXIS = Axis.representatives()['2A']

N_AXES = 20
N_ROUNDS = 5
AXES = []
for i in range(N_AXES):
    AXES.append(STD_AXIS  * MM('r', 5))


start_time = time.time()
for i in range(N_ROUNDS):
    for j in range(N_AXES):
         AXES[j].copy().reduce_G_x0()
end_time = time.time()
run_time = end_time - start_time  # Calculate elapsed time
 
# Optional: Print the average run time
average_time = run_time / (N_AXES * N_ROUNDS) 
print(f"\nAverage run time: {1000 * average_time:.3f} ms")        
    






