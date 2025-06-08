import sys
import os
import time
import shutil
import numpy as np
from argparse import ArgumentParser
import shelve
import time


sys.path.append("axis")
sys.path.append("utilities")
from requirements import check_requirements
check_requirements()

from mat24_orbits import check_recompute, compute_orbits
from mat24_orbits import load_orbits, check_orbits, display_orbits
from suborbits import check_monster_axes
from display_suborbits import display_suborbit_table
from check_all_suborbits import check_suborbits, print_input_for_gap
from check_all_suborbits import call_gap
from cleanup import remove_intermediate_files

CERTIFICATE_PATH = os.path.join("certificates", "axis_certificate.txt")



def parse_args():
    description = ('Display information about G_x0 '
    'orbits of 2A axes. ' 
    )
    usage = "usage: %prog [options]"
    parser = ArgumentParser(usage, description = description)
    parser.add_argument("--all",  dest="all",
        action="store_true",
        help = "Display all information available")
    parser.add_argument("--check-cert",  dest="check_cert",
        action="store_true",
        help = "Check a certificate for computing the number of axes")
    parser.add_argument("--make-cert",  dest="make_cert",
        action="store_true",
        help = "Generate a certificate for computing the number of axes")
    parser.add_argument("--show-Gx0-cent",  dest="show_Gx0_cent",
        action="store_true",
        help = "Display generators for centralizers of G_x0 orbits of axes")
    parser.add_argument("--show-Gx0-orbits",  dest="show_Gx0_orbits",
        action="store_true",
        help = "Display G_x0 orbits of axes and their centralizers")
    parser.add_argument("--show-Gx0-orders",  dest="show_Gx0_orders",
        action="store_true",
        help = "Display orders of centralizers of G_x0 orbits of axes")
    parser.add_argument("--show-N0-orbits",  dest="show_N0_orbits",
        action="store_true",
        help = "Display table of N_0 orbits of axes and their centralizers")
    parser.add_argument("--show-Nx0-orbits",  dest="show_Nx0_orbits",
        action="store_true",
        help = "Display N_x0 orbits of axes, their sizes and G_x0 orbits")
    parser.add_argument("--show-Nxyz-orbits",  dest="show_Nxyz_orbits",
        action="store_true",
        help = "Display Decomposition of N_x0 orbits into N_xyz orbits")
    parser.add_argument("--show-suborbits",  dest="show_suborbits",
        action="store_true",
        help = "Display suborbit diagram for G_x0 orbits of axes")
    parser.add_argument("-r",  dest="recompute", action="store_true",
        help="Recompute all precomputed data")
    parser.add_argument("-t",  dest="latex", action="store_true",
        help = "Display data in format suitable for LaTex (if supported)")
    parser.add_argument("-v",  dest="verbose", action="store_true",
        help="Verbose operation" )
    
    options  = parser.parse_args()
    return options



N_DISPLAY_BLOCKS = 0
def new_block(is_block):
    global N_DISPLAY_BLOCKS
    if N_DISPLAY_BLOCKS:
        print("_" * 78)
    N_DISPLAY_BLOCKS |= bool(is_block)


if __name__ == "__main__":
    options = parse_args()
    d_all = options.all
    if options.recompute:
        remove_intermediate_files()
    recompute = check_recompute(options.recompute)
    if recompute:
        start_time = time.time()
        new_block(recompute)
        d = compute_orbits(n_generators = 10, store = True)
        check_monster_axes()
        from eigenvals_monster import compute_orbits
        compute_orbits()
        from watermark_suborbits import watermark_suborbits
        watermark_suborbits()
        from centralizer_orders import centralizer_orders
        centralizer_orders()
        check_suborbits(check = True, verbose = False)
        t = time.time() - start_time
        T = "Run time for generating tables: %.2f s"
        print(T % t)
        print_input_for_gap()
        gap_ok = call_gap()
    d = load_orbits()   
    check_orbits(d)
    if d_all or options.show_Gx0_cent:
        new_block(all or options.show_Gx0_cent)
        display_orbits(d)
    if d_all or options.show_suborbits:
        new_block(all or options.show_suborbits)
        from display_suborbits import display_suborbit_table
        display_suborbit_table(options.latex)
    if d_all or options.show_Gx0_orbits:
        new_block(all or options.show_G_x0_suborbits)
        from eigenvals_monster import show_eigenvals
        show_eigenvals(options.latex, Nxyz=True)
    if d_all or options.show_Gx0_orders:
        new_block(all or options.show_G_x0_orders)
        from centralizer_orders import centralizer_orders
        centralizer_orders(recompute = False, verbose = 1)

    if d_all or options.show_Nx0_orbits:
        new_block(all or options.show_N_x0_orbits)
        check_suborbits(check = False, verbose = True)
    if d_all or options.show_N0_orbits:
        new_block(all or options.show_N0_orbits)
        from display_N0_suborbits import display_N0_orbits
        display_N0_orbits(options.latex)
    if d_all or options.show_Nxyz_orbits:
        new_block(all or options.show_Nxyz_orbits)
        from display_N0_suborbits import display_Nxyz_orbits
        display_Nxyz_orbits(options.latex)
    if d_all or options.make_cert:
        new_block(all or options.make_cert)
        from make_certificate import make_certificate
        make_certificate(CERTIFICATE_PATH)
    if d_all or options.check_cert:
        new_block(all or options.check_cert)
        from certificates.check_axis_certificate import check_certificate
        check_certificate(CERTIFICATE_PATH)
       


