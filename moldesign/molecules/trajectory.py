from __future__ import print_function, absolute_import, division
from future.builtins import *
from future import standard_library
standard_library.install_aliases()

# Copyright 2017 Autodesk Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import copy
import time

import numpy as np

import moldesign as mdt
from .. import helpers, utils
from .. import units as u
from .molecule import MolecularProperties
from . import toplevel


class Frame(utils.DotDict):
    """ A snapshot of a molecule during its motion. This is really just a dictionary of properties.
    These properties are those accessed as ``molecule.properties``, and can vary
    substantially depending on the origin of the trajectory. They also include relevant dynamical
    data and metadata (such as ``time``, ``momenta``, ``minimization_step``, etc.)

    Properties can be accessed either as attributes (``frame.property_name``) or as keys
    (``frame['property_name']``)

    Some properties with specific meaning are described below:

    Attributes:
        annotation (str): text describing this frame (will be displayed automatically when
            visualized)
        minimization_step (int): for minimization trajectories
        time (u.Scalar[time]): time during a dynamical trajectory

    Example:
        >>> mol = mdt.from_name('benzene')
        >>> mol.set_potential_model(moldesign.models.RHF(basis='3-21g'))
        >>> traj = mol.minimize()
        >>> starting_frame = traj.frames[0]
        >>> assert starting_frame.potential_energy >= traj.frames[-1].potential_energy
        >>> assert starting_frame.minimization_step == 0
    """
    def __init__(self, traj, frameidx):
        self.traj = traj
        self.frameidx = frameidx
        super().__init__(_dynamic=False)
        for key in self.traj.properties:
            self[key] = getattr(traj, key)[self.frameidx]

    def __str__(self):
        return 'Frame %d in trajectory "%s"' % (self.frameidx, self.traj.name)

    def __repr__(self):
        try:
            return '<%s>' % self
        except (KeyError, AttributeError):
            return '<Frame at %x (exception in __repr__)>' % id(self)

    def write(self, *args, **kwargs):
        self.traj._apply_frame(self)
        return self.traj._tempmol.write(*args, **kwargs)

    def as_molecule(self):
        self.traj._apply_frame(self)
        return mdt.Molecule(self.traj._tempmol)


class _TrajAtom(object):
    """ A helper class for querying individual atoms' dynamics
    """
    def __init__(self, traj, index):
        self.traj = traj
        self.index = index
        self.real_atom = self.traj.mol.atoms[self.index]

    @property
    def position(self):
        return self._arrayslice('positions')

    @property
    def momentum(self):
        return self._arrayslice('momenta')

    @property
    def force(self):
        return self._arrayslice('forces')

    def _arrayslice(self, attr):
        return getattr(self.traj, attr)[:, self.index, :]

    def __getattr__(self, item):  # TODO: remove and replace all __getattr__
        if item in ('traj', 'index', 'real_atom'):
            raise AttributeError('_TrajAtom.%s not assigned (pickle issue?)' % item)

        try:  # try to get a time-dependent version
            trajslice = getattr(self.traj, item)
        except AttributeError:  # not time-dependent - look for an atomic property
            return getattr(self.real_atom, item)

        if trajslice[0]['type'] != 'atomic':
            raise ValueError('%s is not an atomic quantity' % item)
        else:
            return u.array([f[self.real_atom] for f in trajslice])

    def __dir__(self):
        attrs = [self.ATOMIC_ARRAYS.get(x, x) for x in dir(self.traj)]
        attrs.extend(dir(self.atom))

    def distance(self, other):
        return mdt.geom.distance(self, other)


class TrajectoryAnalysisMixin(object):
    """
    Organizational class with analysis methods for trajectories
    """
    @property
    def kinetic_energy(self):
        convert_units = True
        energies = []
        for frame in self.frames:
            if 'momenta' in frame:
                energies.append(
                    helpers.kinetic_energy(frame.momenta, self.mol.dim_masses))
            else:
                convert_units = False
                energies.append(None)
        if convert_units:
            arr = u.array(energies)
            return u.default.convert(arr)
        else:
            return energies

    @property
    def kinetic_temperature(self):
        convert_units = True
        temps = []
        energies = self.kinetic_energy
        dof = self.mol.dynamic_dof
        for energy, frame in zip(energies, self.frames):
            if energy is not None:
                temps.append(helpers.kinetic_temperature(energy, dof))
            else:
                convert_units = False
                temps.append(None)
        if convert_units:
            arr = u.array(temps)
            return u.default.convert(arr)
        else:
            return temps

    def distance(self, a1, a2):
        a1, a2 = list(map(self._get_traj_atom, (a1, a2)))
        return mdt.distance(a1, a2)

    def angle(self, a1, a2, a3):
        a1, a2, a3 = list(map(self._get_traj_atom, (a1, a2, a3)))
        return mdt.angle(a1, a2, a3)

    def dihedral(self, a1, a2, a3=None, a4=None):
        a1, a2, a3, a4 = list(map(self._get_traj_atom, (a1, a2, a3, a4)))
        return mdt.dihedral(a1, a2, a3, a4)

    def rmsd(self, atoms=None, reference=None):
        r""" Calculate root-mean-square displacement for each frame in the trajectory.

        The RMSD between times :math:`t` and :math:`t0` is given by

        :math:`\text{RMSD}(t;t_0) =\sqrt{\sum_{i \in \text{atoms}} \left( \mathbf{R}_i(t) -
        \mathbf{R}_i(t_0) \right)^2}`,

        where :math:`\mathbf{R}_i(t)` is the position of atom *i* at time *t*.


        Args:
            atoms (list[mdt.Atom]): list of atoms to calculate the RMSD for (all atoms in the
                ``Molecule``)
            reference (u.Vector[length]): Reference positions for RMSD. (default:
                ``traj.frames[0].positions``)

        Returns:
            u.Vector[length]: list of RMSD displacements for each frame in the trajectory

        """
        if reference is None: refpos = self.frames[0].positions
        else: refpos = reference.positions

        atoms = mdt.utils.if_not_none(atoms, self.mol.atoms)
        indices = np.array([atom.index for atom in atoms])

        rmsds = []
        for f in self.frames:
            diff = (refpos[indices] - f.positions[indices])
            rmsds.append(np.sqrt((diff*diff).sum()/len(atoms)))
        return u.array(rmsds).defunits()


@toplevel
class Trajectory(TrajectoryAnalysisMixin):
    """ A ``Trajectory`` stores information about a molecule's motion and how its properties
    change as it moves.

    A trajectory object contains
      1. a reference to the :class:`moldesign.Molecule` it describes, and
      2. a list of :class:`Frame` objects, each one containing a snapshot of the molecule at a
    particular point in its motion.

    Args:
        mol (moldesign.Molecule): the trajectory will describe the motion of this molecule
        unit_system (u.UnitSystem): convert all attributes to this unit system (default:
            ``moldesign.units.default``)
        first_frame(bool): Create the trajectory's first :class:`Frame` from the molecule's
            current position

    Attributes:
        mol (moldesign.Molecule): the molecule object that this trajectory comes from
        frames (List[Frame]): a list of the trajectory frames in the order they were created
        info (str): text describing this trajectory
        unit_system (u.UnitSystem): convert all attributes to this unit system
    """

    draw = helpers.WidgetMethod('trajectory.draw')
    draw_orbitals = helpers.WidgetMethod('trajectory.draw_orbitals')
    plot = helpers.WidgetMethod('trajectory.plot')

    def __init__(self, mol, unit_system=None, first_frame=False, name=None):
        self._init = True
        self.info = "Trajectory"
        self.frames = []
        self.mol = mol
        self.unit_system = utils.if_not_none(unit_system, mdt.units.default)
        self.properties = utils.DotDict()
        self._reset()
        self._tempmol.dynamic_dof = self.mol.dynamic_dof
        self.name = utils.if_not_none(name, 'untitled')
        if first_frame: self.new_frame()

    def _reset(self):
        self._viz = None
        self._atoms = None
        self._tempmol = mdt.Molecule(self.mol.atoms, copy_atoms=True)


    MOL_ATTRIBUTES = ['positions', 'momenta', 'time']
    """List[str]: Always store these molecular attributes"""

    def copy(self):
        newtraj = copy.copy(self)
        newtraj.frames = self.frames[:]
        newtraj.properties = self.properties.copy()
        newtraj._reset()
        return newtraj

    @property
    def num_frames(self):
        """int: number of frames in this trajectory"""
        return len(self)

    __len__ = mdt.utils.Alias('frames.__len__')
    __iter__ = mdt.utils.Alias('frames.__iter__')

    def __getattr__(self, attr):
        if attr == 'properties' or attr not in self.properties:
            return self.__getattribute__(attr)
        else:
            return self.properties[attr]

    @property
    def atoms(self):
        if self._atoms is None:
            self._atoms = self._make_traj_atoms()
        return self._atoms

    def _make_traj_atoms(self):
        return [_TrajAtom(self, i) for i in range(self.mol.num_atoms)]

    def __str__(self):
        return 'Trajectory for molecule "%s" (%d frames)' % (self.mol, self.num_frames)

    def __repr__(self):
        try:
            return '<%s>' % str(self)
        except (KeyError, AttributeError):
            return '<Trajectory object @ %s (exception in repr)>' % hex(id(self))

    def __add__(self, other):
        newtraj = Trajectory(self.mol, unit_system=self.unit_system)
        for frame in self.frames + other.frames:
            newtraj.new_frame(**frame)
        return newtraj

    def new_frame(self, properties=None, **additional_data):
        """ Create a new frame, EITHER from the parent molecule or from a list of properties

        Args:
            properties (dict): dictionary of properties (i.e. {'positions':[...],
            'potential_energy':...})
            **additional_data (dict): any additional data to be added to the trajectory frame

        Returns:
            int: frame number (0-based)
        """
        # get list of properties for this frame
        props = dict(self.mol.properties)
        if properties is not None:
            props.update(properties)

        for attr in self.MOL_ATTRIBUTES:
            if attr not in props:
                props[attr] = getattr(self.mol, attr)

        props.update(additional_data)

        # add properties to trajectory
        for key, value in props.items():
            if key not in self.properties:
                self._new_property(key, value)
            else:
                proplist = getattr(self, key)
                assert len(proplist) == self.num_frames
                proplist.append(value)

        # backfill missing data with None
        for key in self.properties:
            proplist = self.properties[key]
            if len(proplist) < self.num_frames+1:
                assert len(proplist) == self.num_frames
                try:
                    proplist.append(None)
                except TypeError:
                    newpl = list(proplist)
                    newpl.append(None)
                    self.properties[key] = newpl

        # TODO: less flubby way of keeping track of # of frames
        self.frames.append(Frame(self, self.num_frames))


    def _new_property(self, key, value):
        """ Create a new list of properties for each frame. To facilitate analysis, will try
        to create this list as one of the following classes (in order of preference):
          1) resizeable numpy array
          2) resizeable MdtQuantity array
          3) list

        The list of properties will be backfilled with ``None`` if this property wasn't already
        present
        """
        assert key not in self.properties

        if self.num_frames != 0:
            proplist = [None] * self.num_frames
            proplist.append(value)
        else:
            try:
                proplist = self.unit_system.convert(u.array([value]))
            except TypeError:
                proplist = [value]
            else:
                proplist.make_resizable()
                if proplist.dimensionless:
                    proplist = proplist._magnitude

        self.properties[key] = proplist

    def _get_traj_atom(self, a):
        if a is None:
            return None
        elif isinstance(a, mdt.Atom):
            return self.atoms[a.index]
        else:
            assert isinstance(a, _TrajAtom)

    DONOTAPPLY = set(['kinetic_energy'])

    def _apply_frame(self, frame):
        """
        Reconstruct the underlying molecule with the given frame.
        Right now, any data not passed is ignored, which may result in properties that aren't synced up
        with each other ...
        """
        # TODO: need to prevent multiple things using the _tempmol from conflicting with each other
        m = self._tempmol
        for prop in self.MOL_ATTRIBUTES:
            if (prop not in self.DONOTAPPLY
                and prop in self.properties
                and prop is not None):
                setattr(m, prop, frame[prop])
        m.properties = MolecularProperties(m)
        for attr in frame:
            m.properties[attr] = frame[attr]

    # TODO: need to fix this - it creates import problems
    #@utils.args_from(mdt.converters.write_trajectory,
    #                 allexcept=['traj'],
    #                 append_docstring_description=True)
    def write(self, *args, **kwargs):
        return mdt.fileio.write_trajectory(self, *args, **kwargs)

    def align_orbital_phases(self, reference_frame=None):
        """
        Try to remove orbital sign flips between frames.
        If `reference_frame` is not passed, we'll start with frame 0 and
        align successive pairs of orbitals.
        If `reference_frame` is an int, we'll move forwards and backwards from that frame number.
        Otherwise, we'll try to align every orbital frame to those in reference_frame

        Args:
            reference_frame (int or Frame): ``Frame`` containing the orbitals to align with
                (default: align each frame with the previous one)
        """
        if reference_frame is None:
            iframe = 0
            relative_alignment = True
        elif isinstance(reference_frame, int):
            iframe = reference_frame
            relative_alignment = True
        else:
            relative_alignment = False

        if relative_alignment:
            for i in range(iframe+1, self.num_frames):
                self.frames[i].wfn.align_orbital_phases(
                    self.frames[i-1].wfn)
            for i in range(iframe-1, -1, -1):
                self.frames[i].wfn.align_orbital_phases(
                    self.frames[i+1].wfn)
        else:
            for frame in self.frames:
                frame.wfn.align_orbital_phases(reference_frame.wfn)
