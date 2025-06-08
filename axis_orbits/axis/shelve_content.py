r"""Documentation of the content of the shelve used in this application 

Name of the shelve: mat24_orbit_tables
Subdirectory:       shelve

Entries are documented in the order in that they are created.

Entry name: Lin2Orbits
Created by: mat24_orbits.py
type:       dict: str -> list of pickled instances of class Orbit_Lin2

Function ``load_orbits`` in module ``mat24_orbits`` unpickles these
instances and returns a dictionary 
``str --> instance of class Orbit_Lin2``

Let ``orbit_name`` be the name of a G_x0 orbit of 2A axes. A standard
representative ``axis`` of that orbit is given by 
``mmgroup.axes.Axis.representatives()[orbit_name]``.
The functions in module mat24_orbits.py compute a set of generators
of the centralizer ``C`` of the axis ``axis`` in the group
G_x0. This set of generators is stored in the instance ``y`` of class
``mmgroup.general.Orbit_Lin2``. Here ``y`` is the
value of the dictionary for the key ``orbit_name``.

In ``y`` we store the generators of centralizer ``C`` and also the
operation of ``C`` on the leech lattice mod 2. The methods of class
``Orbit_Lin2`` allow us to retrieve the set of generators of ``C``
and to compute the action of ``C`` on the Leech lattice mod 2.
Object ``y`` also contains representatives of all orbits of type-4
vectors in  the Leech lattice mod 2 under the action of ``C``. 



Entry name: Lin2Samples
Created by: mat24_orbits.py
type:       dict: str -> 2-dimensional numpy array 

Let ``orbit_name`` be the name of a G_x0 orbit of 2A axes. Let ``y``
be the value of "<shelve>["Lin2Orbits"][orbit_name]", as in the
description of entry **Lin2Orbits** of the shelve. Object ``y``
essentially contains the  centralizer ``C`` of the axis given by
``mmgroup.axes.Axis.representatives()[orbit_name]``.
Then method ``y.representatives()`` returns a pair ``(reps, sizes)``,
where ``reps`` is a transversal of the orbits on the type-4 vectors
in the Leech lattice mod 2 under the action of ``C``. If ``reps[i]``
is such a type-4 vector then 
``<shelve>["Lin2Samples"][orbit_name][i]`` is a small array of
type-4 vectors in the same orbit as ``reps[i]``. This array is always
terminated by one or more zero vectors. This array of type-4 vectors
is used for internal tests.


Entry name: mat24_suborbits
Created by: suborbits.py
type:       dict: pair(str, str) -> int

The keys of the dictionary are pairs ``(source, destination)``
indicatings pair of G_x0 orbits. Each G_x0 orbit consists of 8292375
N_x0 orbits, and the trialiy map and its inverse map each member
of an N_x0 orbit to the same pair of destination orbits, if we
consider the latter pair as unordered. Alltogether, there are 
2*8292375 ways how a source G_x0 orbit can be mapped to a destination
orbit. The dictionary value for the key ``(source, destination)``
indicates in how many ways the G_x0 orbit ``source`` is mapped to  
the G_x0 orbit ``destination``. 


Entry name: ORBIT_SIZES
Created by: eigenvals_monster.py
type:       dict: str -> int

Maps the names of the G_x0 orbits of the 2A axes (defined as in entry
"CENTRALIZERS") to their sizes.


Entry name: MAP_SUBORBIT
Created by: watermark_suborbits.py
type:       dict: (watermark) -> int

We number the N_x0 orbits of 2A axes (also called suborbits in file
``watermark_suborbits.py``) with integers. These numbers are ordered
first by the standard order of the G_x0 orbits to which an N_x0 orbit
belongs. In a next step the N_x0 orbits ared orderd by the size of
their centralizers (in decreasing order). Function ``watermark_axis``
in file ``watermark_suborbits.py`` computes some watermark of an axis
that is invariant under the operation of ``N_x0``. This watermark is
a tuple of (tuples of) integers and strings. Dictionary 
``MAP_SUBORBIT`` maps such a watermark of an axis to the number of
the N_x0 orbit of that axis.


Entry name: SUBORBIT_REPRESENTATIVES
Created by: watermark_suborbits.py
type:       list of axes of type Axis

N_x0 orbits of 2A axes are numbered with integers as described in the
previous section. Entry ``i`` of array ``SUBORBIT_REPRESENTATIVES`` 
describes the representative of the N_x0 orbit ``i``. It is a triple
``(orbit_name, entry, v)`` such that the axis representing the
suborbit is ``AXES[orbit_name] * MM0('c', v) ** -1``. Here ``AXES``
is the precomputed dictionary ``mmgroup.axes.Axis.representatives()``
that maps the names of the the name of a :math:`G_{x0}` orbit to 
the axis representing the orbit. 

Let ``lin2_orbits`` be the dictionary obtained from entry 
``Lin2Orbits`` of the shelve. Then 
``reps = lin2_orbits[orbit_name].representatives()[0]`` is the
list of representatives of the orbits of the centralizer of the axis
``AXES[orbit_name]`` on the type-4 vectors in the Leech lattice
(mod2). We have ``v == reps[entry]``.



Entry name: SUBORBIT_SIZES
Created by: watermark_suborbits.py
type:       list of integers

N_x0 orbits of 2A axes are numbered with integers. Entry ``i`` of
array ``SUBORBIT_SIZES`` contains the size of the N_x0 orbit ``i``.



Entry name: SUBORBIT_CENTRALIZERS
Created by: watermark_suborbits.py
type:       list of list of strings

N_x0 orbits of 2A axes are numbered with integers. For each of these 
orbits we compute a representative as in the description of entry
"SUBORBIT_REPRESENTATIVES" of thei shelve. Entry ``i`` of array
``SUBORBIT_CENTRALIZERS`` is a list of strings corresponding
to random elements of the centralizer of the representative of the
the i-th N_x0 orbit. The (random) generators of a centralizer are
modified so that most the first of these generators will be odd,
i.e. in N_x0 \ N_xyz.



Entry name: SUBORBIT_SIZES_2
Created by: check_all_suborbits.py
type:       list of triples ``(e, g, s)``

N_x0 orbits of 2A axes are numbered with integers. Entry ``i`` of
array ``SUBORBIT_SIZES_2`` dscribes the 2 structure of the
centralizer of the i-th N_x0 orbit. Here :math:`N_x0` has structure
:math:`2^[2+11+22}.(M_{24} \times 2)`. The centralizer of the i-th
N_x0 orbit has a corresponding structure :math:`2^e.(G \times 2^s)`.
For the triple  ``(e, g, s)`` the entry ``e`` is 5-tuple describing
the structure of :math:`2^e`, the entry ``g`` is the order of the
group :math:`G`, and entry ``s`` is the exponent of the group 
:math:`2^s`.

Part ``e`` is a tuple of length 5 describing a subgroup of the
2 group :math:`2^{2+11+22}`.  The exponent of the 2 group is
decomposed  as : math:`1+1+11+11+11`, with the i-th component
``e_i`` generated by the following generators:
    
e_0:  :math:`x_{-1}`  
    
e_1:  :math:`x_{\Omega}`  
    
e_2:  :math:`x_{\delta}, \;  \delta \in \mathcal{C}^* \, \mbox{even}`  
    
e_3:  :math:`x_{d}, \;  d \in \mathcal{C}  
    
e_4:  :math:`y_{d}, \;  d \in \mathcal{C}  

For more details, see function ``order_Nx0`` in module ``utilities``.
"""

