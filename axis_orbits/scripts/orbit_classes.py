"""Check representatives of G_x0 orbits of 2A axes

This script checks that the names of the G_x0 orbits of the
2A axes used in the mmgroup package are consistent with the
corresponding names used in [Nor98].
"""


from multiprocessing import Pool, cpu_count
from mmgroup import  MM, MMV
from mmgroup.axes import Axis


def partial_trace(start, end, g, p):
    """Auxiliary function for function ``trace_parallel``."""
    V = MMV(p)
    return sum([int((V([('E', i)]) * g)['E'][i])
          for i in range(start, end)])


def trace_parallel(g, p, n_processes = None):
    """Compute the character of an element of the Monster

    The function returns the character of the element ``g`` of the
    Monster in the 196884-dimensional represention of the Monster
    modulo the value ``p``. The function uses the mmgroup package
    and linear algebra for computing the character. See function
    ``trace_mod_p`` for legal values of the modulus ``p``.

    The optional parameter ``n_processes`` may be used to control 
    the number of parallel processes. 
    """
    if n_processes is None:
        n_processes = max(1, cpu_count() - 1)
    chunk_size = (196884 + n_processes - 1) // n_processes
    jobs = []
    for i in range(0, n_processes):
        i_start = i * chunk_size
        i_end = i_start + chunk_size if i < n_processes - 1 else 196884
        jobs.append((i_start, i_end, g, p))

    with Pool(processes = n_processes) as pool:
        result = pool.starmap(partial_trace, jobs)
    pool.join()
    return sum(result) % p


 
def trace_mod_p(g, p, n_processes = None):
    """Compute order and characters of element of the Monster

    Given an element ``g`` of the Monster, the function returns a
    pair ``(order, characters)``. Here ``order`` is the order of the
    element ``g``. Output  ``characters`` is a dictionary containing
    the character of ``g`` (and, may be, of some powers of ``g``) in
    the representation of the Monster of degree 196883, modulo the
    integer ``p``. More precisely, ``characters[e]`` is the character
    of ``g**e`` (mod ``p``). The value ``characters[1]`` is always
    present in the dictionary. Some more character values may or may
    not be computed. Legal values for modulus ``p`` are
    3, 7, 15, 31, 127, 255. 
    """
    order, chi, _ = MM(g).chi_powers(1,500)
    if chi[1] is not None:
         for e, t in chi.items():
             chi[e] = t % p if t is not None else None
         return order, chi
    else:
        tr = trace_parallel(MM(g), p)
        return (order, {1 : (tr - 1) % p})

 

# The following character information is taken from the ATLAS.
# Here ``CLASSES[cls]`` is a dictionary containing character
# information about the character of degree 196883 for the
# conjugation class ``cl``. If element ``g`` of the Monster is
# of class ``cl`` then ``CLASSES[cl][e]`` is the character
# of ``g**e``. 
CLASSES = {
'2A': {1:4371},
'2B': {1:275},
'4A': {1:275, 2:275},
'4B': {1:51},
'4C': {1:19},
'6A': {1:78},
'6C': {1:14},
'6F': {1:-1},
'8B': {1:11},
'10A': {1:21},
'10B': {1:5},
'12C': {1:6},
}
# The user may easily check that the information provided by
# this dictionary is sufficient to distinguish between the 
# classes mentioned in [Nor98], Table 2.




def check_axes():
    """Check representatives of G_x0 orbits of 2A involutions

    The mmgroup packages contains a precmputed list of 
    representatives of G_x0 orbits of 2A axes. Each of these orbits
    is named by the class of the the element t*x of the Monster
    as in [Nor98], Table 2. Here t is the involution corresponding
    to an axis; and x is the central involution in G_x0.

    Using the character information for the Monster, we check
    that the name of the G_x0 orbits of 2A axes (or involutions)
    in the mmgroup package are compatible with those in [Nor98].
    It suffices to perform the character calculations module 127.
    """
    AXIS_REPS = Axis.representatives()
    x = MM('x', 0x1000)
    p = 127
    for name, ax in AXIS_REPS.items():
        print ("Checking Axis type", name)
        # Let ``t``  be the 2A involution representing the G_x0 orbit
        # called ``name`` in mmgroup.
        t = ax.g_axis
        # Compute order and character of  x * t in he Monster (mod p)
        order, chi = trace_mod_p(x * t, p)
        # Check that order and character information match the
        # corresponding data in [Nor98]. 
        chi_expected = CLASSES[name]
        assert order == int(name[:-1])
        for i in chi_expected:
            assert chi[i] ==  chi_expected[i] % p
        print("Characters are ok")


if __name__ == "__main__":
    check_axes()


   