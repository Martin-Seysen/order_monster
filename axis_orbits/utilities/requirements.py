import sys
import os

######################################################################
# Check if mmgroup and other required packages are installed
######################################################################


def get_mmgroup_version():
    """Return mmgroup version as a triple of ints or None if not present

    """
    try:
        from importlib.metadata import version
        from mmgroup import MM
    except:
        return None
    try:
        return tuple(map(int, version("mmgroup").split(".")))
    except:
        return None


def get_executable():
    """Return name of executable, usually 'python' or 'python3'"""
    return os.path.split(sys.executable)[-1]



def check_mmgroup_version():
    """Check that a suitable version of mmgroup is installed.

    The function raises ModuleNotFoundError and prints a 
    help message if this is not the case.
    """
    version = get_mmgroup_version()
    ok = False
    if version is None:
        print("\nPlease install the mmgroup package!")
        err = "Module mmgroup not found" 
    elif version < (1,0,5):
        print("\nPlease upgrade the mmgroup package to version 1.0.5 or higher!")
        err = "Module mmgroup is outdated" 
    else:
        ok = True
    if not ok:
        py = get_executable()
        print(f"""
For installing or upgrading the mmgroup package, please type in a console:

{py} -m pip install --upgrade mmgroup 

For details we refer to
https://mmgroup.readthedocs.io/en/latest/api.html

""")
        raise ModuleNotFoundError(err)    


def check_other_requirements():
    try:
        import numpy
        import sympy
    except:
        py = get_executable()
        print(f"""
The packages numpy and symmy are required!     

For installing these packages, please type in a console:

{py} -m pip install --upgrade numpy sympy 
""")
        err = "Module numpy or sympy not found"
        raise ModuleNotFoundError(err)    
 
def check_requirements():
    check_mmgroup_version() 
    check_other_requirements()  

