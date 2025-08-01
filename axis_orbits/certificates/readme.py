r"""Accompanying code for 'The Order of the Monster Finite Simple Group'

The main purpose of the paper
'The Order of the Monster Finite Simple Group'
by Gerald Hoehn und Martin Seysen is to give a computational existence
proof of the Monster group M and to compute the order of M.

The most difficult steps in the computation of the order are the proofs
of the correctness of Tables 2 and 4 in the paper. So it should be as
simple as possible to verify the correctness of these tables. We discuss
certificates that can be used for a simple check of the correctness of
these tables.


Counting the axes in the Griess algebra
=======================================

A key ingredient for computing the order of the Monster is
Proposition 3.1 which counts the axes in the Griess algebra.
Table 1 in the paper counts the 2A axis in the twelve different
G_x0-orbits of axis in the Griess algebra. Table 1 can easily be
computed from Table 2 in the paper which contains a transition
matrix for the action of the triality element \tau on the
G_x0-orbits of the axes.

For computing Table 2 we have to find representatives of the
N_x0-orbits on the axes, and elements of M mapping the standard
axis v^+ (defined in the paper) to these representatives. Also we
have to compute the action of the triality element \tau and of its
inverse on these representatives, and to map the images under these
actions to known representative of the orbits.

Finding all these mappings requires sophisticated computations with
the mmgroup package, which are difficult for an independent reviewer
to verify. Once these mappings have been found, they can be recorded
in a file which we will call a certificate. Checking the correctness
of the mappings given by the certificate is much easier than finding
such mappings.

For checking the correctness of a certificate, is suffices if we
can perform the following computations in the Monster and in G_x0:

- Transforming an axis by an element of M or G_x0

- Mapping an element of G_x0 to the automorphism group Co_1
  of the Leech lattice mod 2

- Computing in the natural representation of Co_1 on the
  Leech lattice mod 2

If the correctness of the certificate has been checked then
Table 2 can easily be computed from the data in the certificate.

The script 'check_axis_certificate.py' verifies the correctness
of the certificate 'axis_certificate.txt'; both files are in
subdirectory 'certificates'. The certificate is created by
invoking the python script 'axis.py' with the option
'--make_cert'.


Counting the feasible axes
==========================

Let H be the intersection of the group G_x0 with the centralizer
of the standard axis v^+ as in the paper. Table 3 in the paper counts
the 'feasible axes' in the Griess algebra, where feasible axes are
defined as in the paper. The number of feasible axes must also be
known for computing the order of the Monster. Table 3 can easily be
computed from Table 4 in the paper which contains a transition
matrix for the action of the triality element \tau on the H-orbits
of the feasible axes.

As in the previous section we compute a certificate from which
Table 4 can easily be computed. The requirements for verifying this
certificate are the same as for verifying the certificate in the
previous section.

The script 'check_baby_axis_certificate.py' verifies the
correctness of the certificate 'baby_axis_certificate.txt'
described above; both files are in subdirectory 'certificates'.
The certificate is created by invoking the python script
'baby_axis.py' with the option '--make_cert'.


The word shortening algorithm for the Monster in the mmgroup package
====================================================================

We remark that the word shortening algorithm for the Monster in the
mmgroup package depends in a crucial way on the fact that the list
of G_x0-orbits on axes given by the rows and columns in Table 2 is
complete. The word shortening algorithm for the Monster depends also
on the fact that the list of H-orbits on axes given by Table 4 is
complete.

Verifying the correctness of the certificates in the two previous
yields a proof that the List of G_x0-orbits of axes in given by 
Table 2 and the list of H-orbits of feasible axes in given by Table 4
is complete. As indicated above, the word shortening algorithm is
not required for verifying the correctness of these certificates.

The proof of the correctness of the word shortening algorithm in
[Sey24] refers to tables computed by S.P.Norton and J. Mueller.
Since Table 1 - 4 in our paper contain the same information as
tables computed by S.P.Norton and J. Muller, we obtain a proof of
the correctness of the word shortening algorithm that does not
depend on these tables.


The structure of the certificate
================================

The certificate  'axis_certificate.txt' is a human-readable text file.
Each line in that file is a record that has a structure

tag:  value <g>

``tag`` is a string indicating the meaning of the record.
``value`` is a value which is an int or a string depending on the tag.
``<g>`` is an element of the Monster in mmgroup format.
Some parts of a record may not be present.


We use the following tags:

tag 'axis'  value  <g>

  This record starts the description of an orbit of an axis under G_x0.
  ``value`` is the name of the orbit as in Table 2 in the paper.
  ``g`` an element of M describing the representative ``ax := v^+ * g``
  of the orbit, where ``v^+`` is the standard axis as in the paper.

  The following records contain information about the orbit of the
  axis under G_x0, until a record with tag 'end' occurs.

  In the sequel we let ``ax`` be the axis given by the last recent
  record with tag ``ax``. We let ``value`` be the name of axis ``ax``.


tag 'end'

  This record indicates the end of the description of an axis orbit.


tag 'cent' value <g>

  This record indicates that ``g`` is an element of the centralizer
  ``C(ax)`` of the axis ``ax`` in G_x0.

  If ``value`` is 1 then ``g`` should be considered when computing
  the orbits of ``C(ax)`` on the Leech lattice mod 2. Otherwise this
  is not necessary. For other purposes it may be useful to record
  more generators of ``C(ax)`` in the certificate.

    
tag 'orb' value <g>

  This record describes an N_x0-orbit of axes inside the
  G_x0-orbit of the axis ``ax``. Here ``g`` is an element of G_x0
  such that ``ax1 := ax * g`` is a representative of that orbit.

  A record with tag 'orb' must be followed by record with
  tag 'tau1'.


tag 'tau1' value <g>

  This record must follow a record with tag 'orb'. Let ``ax1``
  be as in the previous record with tag 'orb'. This record
  describes the axis ``ax_t := ax * \tau``, where ``tau``
  is the triality element in M. Then the axis
  ``ax_t * g`` is equal to the axis with name ``value``
  given by one of the records with tag ``axis``.


tag 'tau2' value <g>

  This record must follow a sequence of records with tag 'orb'
  and 'tau1'. Let ``ax1`` be as in the previous record with
  tag 'orb'. This record describes the axis
  ``ax_t := ax * \tau*``, where ``tau*`` is the inverse of
  the triality element in M. Then the axis ``ax_t * g`` is
  equal to the axis with name ``value`` given by one of the
  records with tag ``axis``. 
    


The certificate  'baby_axis_certificate.txt' is a text file of
the same structure as the certificate  'axis_certificate.txt' .
There are just the following differences:

- The certificate describes H-orbits of feasible axes instead of
  G_x0-orbits of axes. Thus an entry ``<g>`` in any record must
  lie in H; i. e. it must centralize the standard axis ``v^+``.

- An entry ``<g>`` in are record with tag 'axis' means that
  the representative of the H-orbit is ``ax := v^- * g``,
  where ``v^-`` is the feasible axis defined in the paper.
"""
