import sys
import os
import time
from collections import defaultdict, OrderedDict
from random import randint, sample
import numpy as np
from argparse import ArgumentParser
import shelve

from mat22_orbits import SHELVE_NAME 
from mat22_orbits import AXES 
ORBITS = list(AXES.keys())



def display_any_table(row_names, col_names, f_values, width, name_width):
    print("%*s" % (name_width+1, ""), end = "")
    for col in col_names:
        print("%*s" % (width, col), end = "")
    print("")
    for row in row_names:
         print("%*s:" % (name_width,row), end = "")
         for col in col_names:
             n = f_values(row, col)
             s = "%*s" % (width, n if n else ".")
             print(s, end = "")
         print("")



def display_any_table_tex(row_names, col_names, f_values, width, name_width):
    print("%*s" % (name_width+1, ""), end = "")
    for col in col_names:
        print(" & ", end = "")
        print("%*s" % (width, col), end = "")
    print(r" \\")
    for row in row_names:
         print("%*s" % (name_width,row), end = "")
         for col in col_names:
             print(" & ", end = "")
             n = f_values(row, col)
             s = "%*s" % (width, n if n else ".")
             print(s, end = "")
         print(r" \\")


HD = r"""Suborbit diagram for H orbits of feasible axes

The table diplays the action of the triality element \tau on the orbits 
H on feasible 2A axes. The entry in row i, column j is the number of
N_xyz orbits of axes in H orbit j that are mapped to H orbit i.
"""


def display_suborbit_table(latex):
    print(HD)
    a = np.zeros((10,12), dtype = np.uint32)
    with shelve.open(SHELVE_NAME) as db:
        d = db["mat22_suborbits"]
    names = ORBITS
    col_tables = (names,)
    f = lambda row, col: d[(col,row)]
    for col_names in col_tables:
        if latex:
            display_any_table_tex(names, col_names, f, 7, 4)
        else:
            display_any_table(names, col_names, f, 7, 4) 
        print("")




def parse_args():
    description = ('Display action of tiality element on representaives of H '
    'orbits of 2A axes. Here H is the intersection of G_x0 and 2.B.' 
    )
    usage = "usage: %prog [options]"
    parser = ArgumentParser(usage, description = description)
    parser.add_argument("-t",  dest="latex", action="store_true",
        help = "Output data in format suitable for LaTex")
    options  = parser.parse_args()
    return options

   

if __name__ == "__main__":
    options = parse_args()
    display_suborbit_table(options.latex)

