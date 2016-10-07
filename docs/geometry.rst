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

    >>> moldesign.angle(atom1, atom2, atom3)
    [...] radians
    >>> moldesign.dihedral(atom1, atom2, atom3, atom4)
    [...] radians

Measuring a dihedral (or twist) angle can be a pain, because they are defined as a function of *four* atomic positions.  For quick measurements, you can pass in only the two central atoms in the dihedral, and MDT will infer the others. Any of the following call signatures will work:

   >>> moldesign.dihedral(a1, a2, a3, a4)
   >>> moldesign.dihedral(a1, a2)
   >>> moldesign.dihedral(bond)


Manipulate geometry
-------------------

WIP



Monitor geometry
----------------

:class:`Geometry monitor <moldesign.geom.GeometryMonitor>` classes can help to track and manipulate of geometric parameters.

   >>> ethylene = mdt.from_smiles('C=C')
   >>> distance_monitor = mdt.DistanceMonitor(ethylene.atoms[0], ethylene.atoms[1])
   >>> distance_monitor.value
   <Quantity(1.51205701815, 'ang')>
   >>> ethylene.draw3d()
.. image:: img/ethyl_init.png
   :scale: 50%

Changing the monitor's value will *change the molecular geometry*:

   >>> distance_monitor.value = 3.0 * u.angstrom
   >>> ethylene.draw()
.. image:: img/ethyl_after.png
   :scale: 50%







Analyze dynamics
----------------

WIP


Constrain geometry
------------------

WIP