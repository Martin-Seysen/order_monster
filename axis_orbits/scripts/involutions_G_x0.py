r"""Show: All involutions in G_{x0} are in class 2A or 2B in the Monster

Here G_{x0} is the maximal subgroup of structure 2^{1+24}.Co_1 of the
Monster. We will show this fact with a brute-force computation uning the
python *mmgroup* package. For the proof we use standard facts about the
Leech lattice (mod 2) and its automorphism group Co_1, and the character
table of Co_1 in the ATLAS [CCN+85], but not any character table
of the Monster or of the group G_{x0}.

An involution in the group G_{x0} of structure 2^{1+24}.Co_1 of the
Monster is either in the subgroup Q_{x0} of structure 2^{1+24} of
G_{x0}, or it maps to one of the classes 2A, 2B, or 2C of involutions
in Co_1. Here classes in Co_1 are denoted as in the ATLAS [CCN+85].

The classes of involutions of G_{x0} contained in Q_{x0} are well
known. They can be described via the the canonical homomorphism from
Q_x0 to the Leech lattice mod 2 given in [Con85], §7. There a vector
in the Leech lattice mod 2 has a type. From Theorem 2 in [Con85] we
conclude that an involution in  $\Qx0$ must map to a vector of type
0, 2, or 4 in the Leech lattice mod 2. G_{x0} is transitive on
elements  of G_{x0} that map to vectors of type 2 or 4, respectivly;
and the unique preimage of a type-0 vector that is an involution in
Q_{x0} is the central involution in  Q_{x0}. The triality element
(obtained as ``MM('t', 1)`` in *mmgroup*) conjugates the central
involution in Q_{x0} to a element of  Q_{x0} corresponding to
a vector of type 4. 

Thus the involutions in Q_{x0} fuse to at most two classes in
the Monster; and a character calculation with *mmgroup* in the
198884-dimensional representation of the Monster shows that these
two classes are different. These two classes of involutions in the
Monster are called 2A and 2B in [CCN+85].

So it suffices to deal with the involutions in G_{x0} that map to
one of the classes 2A, 2B, or 2C in Co_1. It shows that all such 
involutions are in class 2A or 2B of the Monster.

We propose elements I2A, I2B, I2C of G_{x0} that map to the classes
2A, 2B, 2C in Co_1. We check that the squares of these elements are
in Qx0; and we use character calculations in the *mmgroup* package to
show that these elements map to involutions of the appropiate class
in Co_1. Note that the character of degree 299 of Co_1 in the ATLAS
is sufficient to distinguish between elements of Co_1 of classes
1A, 2A, 2B, and 2C. This character is available in *mmgroup*.

Then we use brute force to show that all involutions in the coset
I * Q_{x0}, with I = I2A, I2B, I2C, are in class 2A or 2B in the
Monster. Therefore we use method ``conjugate_involution`` in class
``MM``. This method conjugates any involution in the Monster into
the standard 2A or 2B involution of the Monster with a very high
probability and raises an exception in case of failure. Since this
script should provide a proof, we use this conjugation method as
an oracle, and we recalculate the corresponding conjugation
operations in this script.
"""

from mmgroup import MM, Xsp2_Co1, XLeech2
from mmgroup.bitfunctions import pivot_binary_low, lin_table

ONE = Xsp2_Co1(1)               # neutral element of G_x0
NEG = Xsp2_Co1('x', 0x1000)     # central involution x_{-1} in G_x0

I2A = Xsp2_Co1('y_80fh')        # involution mapping to class 2A in Co_1 
I2C = Xsp2_Co1('y_0ae0h*d_20h') # involution mapping to class 2C in Co_1
I2B = Xsp2_Co1('y_9d0h*d_700h*p_21289010')
                    # I2B squares to x_{-1} and maps to class 2B in Co_1 

def g_centralizer_sign(g):
    """Return centralizer of g in Q_x0 (up to sign)

    Let ``g`` be an element of the group G_x0 given as an
    instance of class ``Xsp2_Co1``. The function returns the
    list of elements of Q_x0 that commute with ``g`` up to
    sign. So this is the list of elements q of Q_x0 such that
    [g, q] is 1 or x_{-1}. Entries in the returned list are
    instances of class ``XLeech2``. The returned list contains
    only one of the elements q, q * x_{-1} of the centralizer. 
    """
    m = g.as_compressed_Co1_bitmatrix()
    m_extended = [int(m[i]) ^ (0x1000001 << i) for i in range(24)] 
    m_map, _ = pivot_binary_low(m_extended)
    m_ker =  [x >> 24 for x in m_map if x & 0xffffff == 0]
    return  [XLeech2(x) for x in lin_table(m_ker)]


def involutions_Q_x0(g):
    """Return list of involutions in the coset g * Q x0

    Here ``g`` must be an element of G_x0 that squares to 1 or x_{-1},
    encoded as an instance of class  ``Xsp2_Co1``. The function
    returns the list of all involutions in the coset  g * Qx0. Entries
    of the list are returned as instances of class ``Xsp2_Co1``.
    """
    involution_list = []
    assert g * g in [ONE, NEG]
    g_square = XLeech2(g * g)
    # Up to sign we have g**2 = 1 and q**2 = 1 for all q in Q_x0.
    # Thus any q in Q_x0 with (g * q)**2 == 1 must commute
    # with g up to sign.
    c_list = g_centralizer_sign(g) # centralizer of g up to sign
    for xl in c_list:
        # Let x = Xsp2_Co1(xl), with xl in class XLeech2.
        # The following condition is equivalent to
        # g**-1 * x * g = x**-1 * g*2, i.e.  (g**-1 * x)**2 == 1
        if xl * g == xl**-1 * g_square:
            xg = Xsp2_Co1(xl) * g
            involution_list += [xg, xg * NEG]
    return involution_list


# Stadard 2A and 2B involutions in the Monster
STD_INVOLUTION = {1: MM('d', [2,3]), 2: MM(NEG)} 

# Static variable for function  reduce_involution():
#    Representatives of classes of involutions in G_x0 found so far
Gx0_REPRESENTATIVES = []

def reduce_involution(g, g_x0_representatives):
    """Conjugate involution to standard 2A or 2B involution

    The function tries to conjugate the involution ``g`` to
    the standard 2A or 2B involution in the Monster. ``g``
    must be an instance of class ``Xsp2_Co1``.
    The function  raises an exception in case of failure.

    Conjugating an involution in G_x0 is much faster than in
    the Monster. So we maintain a list ``g_x0_representatives``
    of representatives of a class of an involution in G_x0.
    We first use method ``conjugate_involution_G_x0`` of class
    ``Xsp2_Co1`` to conjugate ``g1`` to a representative of its
    class in G_x0. If ``g1`` is not in the list
    ``g_x0_representatives`` then we try to conjugate ``g1`` to
    the standard 2A or 2B involution in the Monster; and we
    append ``g1`` to that list. This final step can be dropped
    if ``g1`` is already in that list.
    """
    global Gx0_REPRESENTATIVES
    assert g*g == ONE
    # We use the methods for conjugating involutions in classes
    # Xsp2_Co1 and MM as oracles for finding elements conjugating
    # ``g`` to the representatve of its class. But in this function
    # we check that these conjugations work as expected. We first 
    # conjugate in G_x0, since this is much faster than in M.  
    _, h = g.conjugate_involution_G_x0()
    g1 = g**h
    if g1 not in g_x0_representatives:
        # Then conjugate g1 in M to the representative of class
        # 2A or 2B in the Monster M.
        mg1 = MM(g1)
        i, h = mg1.conjugate_involution()
        # check that h conjugates g to a standard involution in M
        assert mg1**h == STD_INVOLUTION[i]
        # finally, append g1 to the list ``g_x0_representatives``
        g_x0_representatives.append(g1)  

    
def check_involution(name, g):
    """Check all involutions in the coset g * Q_x0

    Here ``g`` must be an element of G_x0 that squares to
    1 or x_{-1}. The function checks that all involutions in
    the coset g * Q_x0 are 2A or 2B involutions in the Monster.
    It raises an exception if this is not the cae.
    """
    # Store all involutions in g * Q_x0 in ``involution_list``
    involution_list = involutions_Q_x0(g)
    print("A preimage of Co_1 class %s in G_x0 contains %d involutions"  
         % (name, len(involution_list)))
    if len(involution_list) > 0:
        # Check that all involutions in the list are in class 2A or 2B
        # in the Monster. We also compute a list ``representatives``
        # such that every involution in ``involution_list`` is 
        # conjugate in G_x0 to an  involution in ``representatives``.
        representatives = []
        for inv in involution_list:
            reduce_involution(inv, representatives)
        print("All involutions are in class 2A or 2B in the Monster")
        # Next show that the involutions in the list ``representatives``
        # can be distiguished by their characters available in mmgroup.
        # Hence these involutions are in different classes in G_x0. 
        characters = []
        G_x0_classes = {1:0, 2:0}
        for g in representatives:
            chi = g.chi_G_x0()
            assert chi not in characters
            characters.append(chi)
            # Count G_x0 classes fusing to class 2A and 2B in the 
            # Monster in values G_x0_classes[1] and G_x0_classes[2].  
            G_x0_classes[g.conjugate_involution()[0]] += 1
        print("""The preimage in G_x0 has %d classes fusing to class 2A
and %d classes fusing to class 2B in the Monster"""
            % (G_x0_classes[1], G_x0_classes[2]))
    print("") 



# The character of degree 299 of Co_1 applied to the involutions
# in Co_1, taken form the ATLAS. 
CHI_299 = {'2A' : 43, '2B' : -13, '2C' : 11}

if __name__ == "__main__":
    for name, inv in zip(['2A', '2B', '2C'], [I2A, I2B, I2C]):
        assert inv * inv in [ONE, NEG]
        assert inv.chi_G_x0()[1] == CHI_299[name], inv.chi_G_x0()
        check_involution(name, inv)

