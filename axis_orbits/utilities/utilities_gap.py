import sys
import os
from collections import defaultdict, OrderedDict
import re
import numpy as np
import subprocess

from mmgroup import Xsp2_Co1, AutPL
from mmgroup.general import Orbit_Lin2
from mmgroup.mat24 import bw24, gcode_to_vect, cocode_weight
from mmgroup.mat24 import vect_to_list, spread_b24, syndrome, lsbit24
from mmgroup.mat24 import cocode_syndrome
from mmgroup.bitfunctions import iter_bitweight


from utilities import MM_to_GAP, is_Nx0_odd, sub_structure_description

######################################################################
# Write GAP code for computing the 2 structure of a subgroup
######################################################################


def sample_upper_central_series(f):
    """Sample for upper central series of a 2 group.

    The function writes GAP code to demonstrate the computation of
    an upper central series of a well-known 2 group to file ``f``.
    """
    print(
"""
Q8 := SmallGroup(8,4);;
UCS := UpperCentralSeriesOfGroup(Q8);;
Print("Example: The quaternion group Q8 has upper central series\\n");
for i in [1..Length(UCS)-1] do
  Print(StructureDescription(FactorGroup(UCS[i],UCS[i+1])), "\\n");
od;;
Print("and order ", Order(Q8), "\\n");
""", file = f)




def subgroups_o2(name,  f):
    """

    For checking if a group extension plits, see:
    https://math.stackexchange.com/questions/3991118/is-there-a-way-to-check-if-a-group-extension-is-split-or-not-using-gap

    For function ComplementClassesRepresentatives,
    see GAP documentation Ch. 39.11-6:
    https://docs.gap-system.org/doc/ref/chap39.html#X804F0F037F06E25E
    """
    print(
f"""O2{name} := PCore({name}, 2);;
UCS := UpperCentralSeriesOfGroup(O2{name});;
Print("O2{name} has upper central series:\\n");
for i in [1..Length(UCS)-1] do
  Print(StructureDescription(FactorGroup(UCS[i],UCS[i+1])), "\\n");
od;;
Print("and order ", Order(O2{name}), "\\n");
Print("{name}/O2{name} = ", StructureDescription(FactorGroup({name}, O2{name}:nice)), "\\n");
SPLIT := Length(ComplementClassesRepresentatives({name},O2{name})) > 0;;
Print("split = ", SPLIT, "\\n");
""", file = f)





######################################################################
# Create input for GAP for computing order and 2-structure of a group
######################################################################


def create_input_for_gap(centralizers, orders, f_in):
    with open(f_in,"wt") as f:
        print("# This is a GAP program.")
        print('LoadPackage("smallgrp");', file = f);
        sample_upper_central_series(f)
        #std_centralizers_gap(f)
        for i, (c, order) in enumerate(zip(centralizers, orders)):
            small_order = order <= 1000 and order % 256 != 0
            small_order_nice = order <= 1000 and order % 256 != 0
            small_order_nice |= order <= 10000 and order % 128 != 0
            #tag, lst = find_substructure(c) if sub else ("", [])
            #conj = AutPL(0, zip(lst, range(5)), 0) if tag else None
            odd = is_Nx0_odd(c[0])
            g = MM_to_GAP(c[odd:])
            #print(i, order, small_order_nice)
            nice = ":nice" if small_order_nice  else ""
            print(f"""G{i} := Group(
{g});;
# G{i} should have order {order}
Print("G{i} = ", StructureDescription(G{i}{nice}), "\\n");
Print("order = ", Order(G{i}), "\\n");
Print("small = {str(bool(small_order_nice)).lower()}\\n"); 
""", file = f, end = "")
            subgroups_o2("G"+str(i),  f)
            if small_order:
                print(f'Print("Id = ", IdSmallGroup(G{i}), "\\n");',
                    file = f)
            print("", file = f)


######################################################################
# Parse output lines from GAP
######################################################################

def f_o2c(slist):
    data =  [y.strip() for y in re.split('x', slist[0])]
    d = defaultdict(int)
    for y in data:
        d[int(y[1:])] += 1
    return [ [(i, d[i])  for i in sorted(d)] ]


MATCHES = [
# match "G<n> = <structure_sedcription>"
[re.compile(r"G(\d+)\s*=\s*(.+)$"), "g", "ds"], 
# match "order = <order>"
[re.compile(r"order\s+=\s+(\d+)"), "order", "d"],
# match "Id = <Id>"
[re.compile(r"Id\s*=\s*(.+)$") , "id", "s"],
# match "small = <true|false>"
[re.compile(r"small\s*=\s*(true|false)\s*$"), "small", "b"],
# match "split = <true|false>"
[re.compile(r"split\s*=\s*(true|false)\s*$"), "split", "b"],
# match "O2G<n>\s+has\s+upper\s+central"
[re.compile(r"O2G(\d+)\s+has\s+upper\s+central"), "o2", "d"],
# match "C<n> [x C<n>]"
[re.compile(r"\s*(C\d+(\s*x\s*C\d+)*)"), "o2c", f_o2c],
# match "and\s+\order\s+(\d+)"
[re.compile(r"\s*and\s+order\s+(\d+)"), "o2order", "d"],
# match "G<n>/O2G<n> = <data>"
[re.compile(r"G(\d+)/O2G(\d+)\s*=\s*(.+)$"), "o2f", "dds"],
]



def iter_parse_gap_output(filename, verbose = 1):
    for s in open(filename, "rt").readlines():
        if verbose:
            print(s, end = "")
        for m, tag, fmt in MATCHES:
            mm = m.match(s)
            if mm:
               data = [tag]
               if isinstance(fmt, str):
                   for fc, x in zip(fmt, mm.groups()):
                      if fc == "b":
                          data.append(eval(x[0].upper() + x[1:]))
                      elif fc == "d": 
                          data.append(int(x))
                      elif fc == "s": 
                          data.append(x.strip())
               else:
                   data += fmt(mm.groups())
               yield data
               break


OP_TYPE_DICT = {
  '.': 1, ':' : 2, 'x' : 3,
}          



######################################################################
# Convert GAP output for the 2-structure of a group to a string
######################################################################



def str_o2_factor(s, short):
    def as_exp(b, e):
        if e == 1:
            return str(b)
        elif e < 10:
            return "%s^%s" % (b, e)
        else:
            return "%s^{%s}"  % (b, e)
    short = short and len(s) == 1 and s[0][0] == 2
    if short:
        return str(s[0][1]), short
    elif len(s) == 1:
        return as_exp(*s[0]), False
    else:
        data = [as_exp(*x) for x in s]
        return "(" + r" \times ".join(data) + ")", False


def str_o2_subgroup(c_list):
    if len(c_list) == 0:
        return "1"
    short = True
    shortlist, longlist = [],[]
    for s in c_list:
        s1, short = str_o2_factor(s, short)
        if short:
            shortlist.append(s1)
        else:
            longlist.append(s1)
    s_long = ".".join(longlist)
    #print("2222<%s><%s>" % (shortlist, s_long))
    if len(shortlist) == 0:
        return s_long
    if len(shortlist) == 1 and shortlist[0] == '1':
        s_short = "2"
    else:
        s_short = "2^{%s}" % "+".join(shortlist)
    if len(longlist) == 0:
        return s_short
    return s_short + "." + s_long
   



######################################################################
# Store output of GAP in an object of class GapInfo
######################################################################




def structure_op_type(s):
    op_type = 'a'
    bra = 0
    precedence = 0
    for i, c in enumerate(s):
        bra += c == '('
        if bra == 0 and c == ' ' and s[i+1:i+2] in OP_TYPE_DICT: 
            op = s[i+1:i+2]
            prec = OP_TYPE_DICT[op]
            if prec > precedence:        
                op_type = op
                precedence = prec
        bra -= c == ')'
    return op_type       


def str_complexity(s):
    if s is None:
        return (1000,1000,1000,1000)
    bra = op = op1 = exp = 0
    for c in s:
        bra += c == '('
        op += c in ".:x"
        op1 += c in ".:"
        exp += c == '{'
    return (bra, op, op1, exp, len(s))



class GapInfo:
    no = None         # current number of group G computed in GAP
    order = None      # order of that group G
    id = None         # GAP id of G (if G is small)
    small = None      # True if G s small enough for a nice structure description
    structure = None  # structure description of G
    o2 = None         # upper central series of O_2(G)
    o2_str = None     # structure of O_2(G) as string
    o2_factor = None  # structure description of G/O_2(G)
    split = None      # True if extension of O_2(G) by G splits

    def o2_structure_description(self):
        if self.o2_str is None or self.o2_factor is None:
            return None
        if self.o2_str == "1":
            return self.o2_factor
        if self.o2_factor == "1":
            return self.o2_str
        factor_type = structure_op_type(self.o2_factor)
        if self.split:
            bra = factor_type in "x"
        else:
            bra = factor_type in "x:."
        factor = "(" + self.o2_factor + ")" if bra else self.o2_factor
        op = " : " if self.split else " . "
        return self.o2_str + op + factor
        

    def structure_description(self):
        s_o2 = self.o2_structure_description()
        if self.small:
            s = self.structure
            if str_complexity(s) <= str_complexity(s_o2):
                return sub_structure_description(s)
            else:
                return sub_structure_description(s_o2)
        elif s_o2 is not None:
            return sub_structure_description(s_o2)
        return sub_structure_description(s)   
            
    


######################################################################
# Store GAP output in an array of objects of class GapInfo
######################################################################




def parse_gap_output(filename, verbose = 0):
    data = []
    no = None
    c_list = o2_prod = None
    for s in iter_parse_gap_output(filename, verbose > 1):
        if verbose > 1:
            print(s)
        tag = s[0]
        if tag == 'g':
            no = len(data)
            g = GapInfo()
            data.append(g)
            g.no, g.structure = s[1:]
            if verbose:
                print("\nstr(G%d) =" % no, g.structure)
            assert g.no == no, (g.no, no)
        if no is None:
            continue
        if tag in ['order', 'small', 'split']:
            setattr(g, tag, s[1])
            assert g.no == no, (g.no, no)
            if verbose:
                print("%s(G%d) =" % (tag,no), s[1])
        elif tag == 'o2f':
            assert s[1] == s[2] == no
            g.o2_factor = s[3]
            if verbose:
                print("F_2(G%d) = %s" % (no, g.o2_factor))
        elif tag == 'o2':
            assert s[1] == no
            c_list, o2_prod = [], 1
        elif tag == 'o2c':
            c_list.append(s[1])
            if verbose > 1:
                print(s[1])
            for (x, e)  in s[1]:
                o2_prod *= x ** e
        elif tag == 'o2order':
            assert s[1] == o2_prod, (s[1] , o2_prod)
            rev =  list(reversed(c_list))
            g.o2 = rev
            g.o2_str = str_o2_subgroup(rev)
            if verbose:
                print("O_2(G%d) =" % no, g.o2_str)
            c_list = o2_prod = None
    if verbose:
        print("")
    d = OrderedDict()
    for x in data:
        d[x.no] = x
    return d          
           
  

######################################################################
# Run a GAP program
######################################################################




def run_process(command, args=None, dir=None, input_file=None, output_file=None):
    # Ensure args is a list if provided
    if args is None:
        args = []
    
    # Check if running on Windows
    is_windows = sys.platform.startswith("win")
    
    # Open input and output files if specified
    if dir:
        input_file = os.path.join(dir, input_file)
        output_file = os.path.join(dir, output_file)
    input_handle = open(input_file, 'r') if input_file else None
    output_handle = open(output_file, 'w') if output_file else None
    print("Running gap <%s >%s" % (input_file, output_file))
    try:
        # Run the process with redirection
        process = subprocess.run([command] + args, stdin=input_handle,
            stdout=output_handle, text=True, shell=is_windows)
        
        # Optionally, print the output if not redirected to a file
        if output_file is None:
            print(process)
    finally:
        # Close file handles if they were opened
        if input_handle:
            input_handle.close()
        if output_handle:
            output_handle.close()




def run_gap(input_file=None, output_file=None, dir = None):
    try:
        prog = f"gap -g <{input_file} >{output_file}"
        if sys.platform.startswith("win"):
            prog = f"wsl gap -g <{input_file} >{output_file}"
            run_process("wsl", ["gap", "-b", "-q"], dir, input_file, output_file)
        else:
            run_process("gap", ["-b", "-q"], dir, input_file, output_file)
        print("gap terminated successfully")
        return True
    except:
        print(f"""
Launching GAP failed!
Please switch to directory "{dir}" and run
{prog}
""")
        raise
        return False

    



######################################################################
# Test function ``parse_gap_output`` of this module
######################################################################


if __name__ == "__main__":
    d = parse_gap_output('Nx0_orbit_structure.txt', verbose = 1)
    for x in d.values():
        print("G%-3d = %-31s = %-30s" % (x.no,
           x.structure_description(), x.structure  ))
    #1/0


