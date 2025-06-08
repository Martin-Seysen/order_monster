r"""Compute eigenspaces of an axis in the Griess algebra

Yet to be documented
"""

from collections import defaultdict
from mmgroup import  XLeech2



def count_scalar_products_short():
    """Count some scalar products of short vectors in the Leech lattice

    Here we write L for the Leech lattice and L2 for the Leech lattice
    mod 2. A shortest vector in L2 has a unique shortest preimage in L
    up to sign. Thus the scalar product <v, w> of two shortest vectors
    v, w in L2 may be defined as the absolute value of the scalar
    product of their preimages in L; this may be 0, 1, 2, or 4. If v
    or w is not short in L2 then the scalar product <v, w> is defined
    modulo 2 only.
    
    Let \Omega be the standard frame in L2; and let \beta by a vector
    x_{ij} in the notation of [Con85]. Then \beta can be characterized
    by the property that \beta + \Omega is short.

    We want to count the scalar products <v, \Omega> (mod 2) and
    <v, \beta>  for all short vectors v in L2. The function returns
    a dictionary d. The key of d are pairs of integers (s, m).
    Then d[(s,m)] is the number of short vectors v in L2 with
    <v, \Omega> = m (mod 2) and <v, \beta> = s.
    """
    # mmgroup supports computation in the double cover Q_x0 of L2 only. 
    OMEGA = XLeech2('Omega') # Preimage in Q_x0 of standard frame in L2
    BETA = XLeech2(0, [2,3]) # Preimage in Q_x0 of short vector x_{2,3}
    d = defaultdict(int)     # Result of the computation
    # Let the type of a vector in L or L2 be as in [CS99], Ch. 10.3.3.
    assert BETA.type == (BETA * OMEGA).type == 2
    assert OMEGA.type == 4
    NUMBER_OF_SHORT_VECTORS = 98280
    for i in range(NUMBER_OF_SHORT_VECTORS):
        # Let v be (a preimage in Q_x0 of) the i-th short vector in L2.
        v = XLeech2('short', i)
        # Compute scalar product <v, BETA> as
        # type(v) + type(BETA) - type(v * BETA).
        skalprod_beta = 4 - (v * BETA).type
        assert skalprod_beta in [0, 1, 2, 4]
        # Compute the scalar product <v, Omega> modulo 2.
        skalprod_Omega = (v * OMEGA).type % 2
        # Count the short vector v with these scalar products
        d[skalprod_beta, skalprod_Omega] += 1
    assert d[(4,0)] == 1  # This counts the short vector BETA only 
    assert d[(4,1)] == 0  # Scalar product <BETA, OMEGA> is even
    # There must be a pairing between v and v + BETA
    # for all v in L2 with <v, BETA> = 2.
    assert d[(2,0)] % 2 == d[(2,1)] % 2 == 0
    return dict(d)


d = count_scalar_products_short()

print("""Scalar products of vectors v in Lambda / 2 Lambda:
               <v,Omega>
             even     odd
<v,beta>""")
for prod in [0, 1, 2, 4]:
    print("   %2d       %5d   %5d" % (prod, d[(prod,0)], d[(prod,1)]))
print()

DIM_A = {'n_0': 1 + 24 * 23 // 2, 'n_2': 23}

n_16 = 1
n_0 = DIM_A['n_0'] + d[0,0] + 3 * d[0,1] + d[2,0] // 2 + 3 * d[2,1] // 2
n_2 = DIM_A['n_2'] + d[2,0] // 2 + 3 * d[2,1] // 2
n_half =  d[1,0] + 3 * d[1,1]
n_sum = n_16 + n_0 + n_2 + n_half
assert n_sum == 196884

print("\nDimensions of Eigenspaces of ad(axis) in the Griess algebra")
print("Eigenvalue:     16      0      2    1/2")
print("Dimension:", ("%7d"*4) % (n_16, n_0, n_2, n_half))  
  
    
