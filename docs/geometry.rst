Geometry tools
==============


Measure geometry
----------------

The :mod:`moldesign.geom` module contains a variety of methods for measuring (and
manipulate geometry).

You can get the distance between any two atoms with the
:meth:`atom.distance <moldesign.Atom.distance>` method.

  >>> atom1.distance(atom2)
  [...] angstrom

Bond angles and dihedral (twist) angles can be measured using the :meth:`moldesign.geom.angle`
and :meth:`moldesign.geom.dihedral` methods:

    >>> moldesign.geom.angle(atom1, atom2, atom3)
    [...] radians
    >>> moldesign.geom.dihedral(atom1, atom2, atom3, atom4)
    [...] radians


Manipulate geometry
-------------------

WIP


Analyze dynamics
----------------

WIP


Constrain geometry
------------------

WIP