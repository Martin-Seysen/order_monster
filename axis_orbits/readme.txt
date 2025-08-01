Introduction
============
This is the accompanying Python code for the paper
"The Order of the Monster Finite Simple Group"
by Gerald Hoehn and Martin Seysen.

This accompanying code is available from:
https://github.com/Martin-Seysen/order_monster

Instructions for installation and test are given on that web page.


Running the scripts axis.py and baby_axis.py
============================================

The scripts axis.py and baby_axis.py compute the Tables 1 - 6
in the paper. Here axis.py deals with the axes, and baby_axis.py
deals with the feasible axes in the Griess Algebra, as described
in the paper.

Each table in the paper can be generated with one of these
scripts and one of the following options:

Option                    Table computed by script
                          axis.py      baby_axis.py

--show-Gx0-orbits         Table 1      Table 3
--show-suborbits          Table 2      Table 4
--show-N0-orbits          Table 5      Table 6

So e.g. Table 2 can be computed with:

python3 axes.py --show-suborbits

In a Windows system the string "python3" should be replaced by "python".


Apart from these options, the following options are available:

 -h, --help          Show a help message containing all possible options
 --show-Gx0-orders   Display orders of centralizers of G_x0-orbits on axes
 --show-Nx0-orbits   Display N_x0-orbits on axes, their sizes and G_x0-orbits
 --show-Nxyz-orbits  Display Decomposition of N_x0-orbits into N_xyz-orbits
 -t                  Display data in format suitable for LaTex (if supported)
 -r                  When the script is called for the first time it computes
                     large internal tables to speed up subsequent calls.
                     This option forces a recomputation of these tables.


Cleaning up
===========

For cleaning up, the script 'cleanup.py' can be run with the following
options:

 -r          remove all intermediate files and internal tables
 -z          zip the source files to the file 'axis_orbits.zip'



Some auxiliary scripts
======================

The following scripts perform auxiliary tasks related to the
information shown in the paper.

For running these scripts, switch to subdirectory 'scripts',
and run any of the following scripts.


eigenspaces_axis.py:

   This script computes the dimensions of the eigenspaces of ad(v)
   for a fixed axis v, as required for the table in the proof of
   Lemma 2.4. Therefore it runs through the standard basis of
   the subspace 98280_x of the Griess Algebra, and computes the
   required information for each of the basis vectors.


involutions_G_x0.py:

   The script enumerates the classes of involutions in the group
   G_x0 using the mmgroup package and brute force.
   This enumeration is mentioned in Appendix D.


orbit_classes.py:

   In Table 1 in the paper the G_x0-orbits on the axes are
   labeled by the name of the class of an element t * x in the
   Monster, where x is the central Involution in G_x0, and t is
   a 2A involution corresponding to any axis in that orbit.
   These names are given in ATLAS notation.

   A transversal V of the G_x0-orbits on the axes is computed
   with mmgroup. Let A be the list of the elements t * x, where
   x runs through the involutions corresponding to the axes in V.
   The script computes the characters of the entries of A in the
   Griess Algebra with mmgroup.

   These characters are used to identify the classes of t * x
   in the Monster in ATLAS Notation.

   The script raises an exception if any of the twelve names in
   Table 1 is not correct.


Certificates
============

The main purpose of the paper is the computation of the order
of the Monster. Here the most difficult step is the computation
of Table 2 in the paper. If this table has been computed then the
index |M : 2.B| which is equal to the number of the axes in the
Monster M can easily be computed.

For computing Table 2, some rather difficult computations with
the mmgroup package must be done. On the other hand, checking
the correctness of Table 2 should be as simple as possible.
Therefore it makes sense to store auxiliary information in
a file in human-readable form that can be used to simplify the
proof of the correctness of the table. We call such a file
a 'certificate'. The description of the structure of the
certificate and a simple program for proving the correctness
of Table 2 is contained in subdirectory 'certificates'.

To generate a certificate one may call

python3 axis.py --make-cert

in the main directory. Then a certificate is generated and stored
in subdirectory 'certificate'. To check the correctness of the
certificate one may call:

python3 axis.py --check-cert

More details are given in file 'certificates/readme.py'.

Thus for an independent verification of the correctness of
Table 2 it suffices to study the certificates and algorithms
given in subdirectory 'certificates' of the accompanying code.

Regarding computations in the Monster group, it suffices if
the verifier has the capability to transform a vector in
the Griess Algebra with an element of the Monster.

A similar process is implemented for computing Table 4.
Here a certificate for containing auxiliary information
for computing Table 4 can also be generated.

To generate such a certificate one may call

python3 baby_axis.py --make-cert

To check the correctness of the certificate one may call:

python3 baby_axis.py --check-cert
 
For checking the correctness of Tables 2 and 4 in the paper it
suffices to generate the certificates mentioned in this section
and to check and execute the scripts
'check_axis_certificate.py' and 'check_baby_axis_certificate.py'
in subdirectory 'certificates'.
