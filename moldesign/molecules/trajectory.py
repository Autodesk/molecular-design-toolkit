# Copyright 2016 Autodesk Inc.
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
import os
import time
from cStringIO import StringIO

import numpy as np

import moldesign as mdt
from moldesign import helpers, utils
from moldesign import units as u

from .molecule import MolecularProperties
from . import toplevel


class SubSelection(object):
    """
    Descriptors to get bits of the trajectory
    trajectory.atoms[3].position -> array of positions
    trajectory.atoms[5].distance( mol.chains['B'].residues[5] ) -> array of distances
    trajectory.chains['A'].residue[23].com -> array of COMs
    NOT IMPLEMENTED YET
    """


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
        >>> mol.set_potential_model(moldesign.methods.models.RHF(basis='3-21g'))
        >>> traj = mol.minimize()
        >>> starting_frame = traj.frames[0]
        >>> assert starting_frame.potential_energy >= traj.frames[-1].potential_energy
        >>> assert starting_frame.minimization_step == 0
    """
    pass


class _TrajAtom(object):
    """ A helper class for querying individual atoms' dynamics
    """

    ATOMIC_ARRAYS = {'position': 'positions',
                     'momentum': 'momenta',
                     'force': 'forces'}

    def __init__(self, traj, index):
        self.traj = traj
        self.index = index
        self.real_atom = self.traj.mol.atoms[self.index]

    def __getattr__(self, item):
        if item in ('traj', 'index', 'real_atom'):
            raise AttributeError('_TrajAtom.%s not assigned (pickle issue?)' % item)

        if item in self.ATOMIC_ARRAYS:
            is_array = True
            item = self.ATOMIC_ARRAYS[item]
        else:
            is_array = False

        try:  # try to get a time-dependent version
            trajslice = getattr(self.traj, item)
        except AttributeError:  # not time-dependent - look for an atomic property
            return getattr(self.real_atom, item)

        if is_array:
            return trajslice[:, self.index, :]
        else:
            if trajslice[0]['type'] != 'atomic':
                raise ValueError('%s is not an atomic quantity' % item)
            else:
                return u.array([f[self.real_atom] for f in trajslice])

    def __dir__(self):
        attrs = [self.ATOMIC_ARRAYS.get(x, x) for x in dir(self.traj)]
        attrs.extend(dir(self.atom))

    def distance(self, other):
        return mdt.geom.distance(self, other)


@toplevel
class Trajectory(object):
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
    def __init__(self, mol, unit_system=None, first_frame=False):
        self._init = True
        self.info = "Trajectory"
        self.frames = []
        self.mol = mol
        self.unit_system = utils.if_not_none(unit_system, mdt.units.default)
        self._property_keys = None
        self._tempmol = mdt.Molecule(self.mol.atoms, copy_atoms=True)
        self._tempmol.dynamic_dof = self.mol.dynamic_dof
        self._viz = None
        self.atoms = [_TrajAtom(self, i) for i in xrange(self.mol.num_atoms)]
        if first_frame: self.new_frame()

    # TODO: the current implementation does not support cases where the molecular topology changes
    # TODO: allow caching to disk for very large trajectories

    MOL_ATTRIBUTES = ['positions', 'momenta', 'time']
    """List[str]: Always store these molecular attributes"""

    def __getstate__(self):
        return self.__dict__.copy()

    def __setstate__(self, d):
        self.__dict__.update(d)

    @property
    def num_frames(self):
        """int: number of frames in this trajectory"""
        return len(self)

    def __len__(self):
        """overrides len(trajectory) to return number of frames"""
        return len(self.frames)

    @utils.kwargs_from(mdt.widgets.trajectory.TrajectoryViewer)
    def draw3d(self, **kwargs):
        """TrajectoryViewer: create a trajectory visualization

        Args:
            **kwargs (dict): kwargs for :class:`moldesign.widgets.trajectory.TrajectoryViewer`
        """
        self._viz = mdt.widgets.trajectory.TrajectoryViewer(self, **kwargs)
        return self._viz
    draw = draw3d  # synonym for backwards compatibility

    def draw_orbitals(self, align=True):
        """TrajectoryOrbViewer: create a trajectory visualization"""
        from moldesign import widgets
        for frame in self:
            if 'wfn' not in frame:
                raise ValueError("Can't draw orbitals - orbital information missing in at least "
                                 "one frame. It must be calculated with a QM method.")

        if align: self.align_orbital_phases()
        self._viz = widgets.trajectory.TrajectoryOrbViewer(self)
        return self._viz

    def __str__(self):
        return 'Trajectory for molecule "%s" (%d frames)' % (self.mol, self.num_frames)

    def __repr__(self):
        try:
            return '<%s>' % str(self)
        except Exception:
            return '<Trajectory object @ %s (exception in repr)>' % hex(id(self))

    def __add__(self, other):
        newtraj = Trajectory(self.mol, unit_system=self.unit_system)
        newtraj.frames = self.frames + other.frames
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
        # TODO: callbacks to update a status display - allows monitoring a running simulation
        if properties is None:
            new_frame = Frame()
            for attr in self.MOL_ATTRIBUTES:
                val = getattr(self.mol, attr)
                try:  # take numpy arrays' values, not a reference
                    val = val.copy()
                except AttributeError:
                    pass
                if val is not None: new_frame[attr] = val
            new_frame.update(self.mol.properties)
        else:
            new_frame = Frame(properties)
        for key, val in additional_data.iteritems():
            assert key not in new_frame, (
                "Can't overwrite molecule's properties with additional_data key %s" % key)
            new_frame[key] = val

        self._update_property_keys(new_frame)
        self.frames.append(new_frame)


    def _get_traj_atom(self, a):
        if a is None:
            return None
        elif isinstance(a, mdt.Atom):
            return self.atoms[a.index]
        else:
            assert isinstance(a, _TrajAtom)

    def distance(self, a1, a2):
        a1, a2 = map(self._get_traj_atom, (a1, a2))
        return mdt.distance(a1, a2)

    def angle(self, a1, a2, a3):
        a1, a2, a3 = map(self._get_traj_atom, (a1, a2, a3))
        return mdt.angle(a1, a2, a3)

    def dihedral(self, a1, a2, a3=None, a4=None):
        a1, a2, a3, a4 = map(self._get_traj_atom, (a1, a2, a3, a4))
        return mdt.dihedral(a1, a2, a3, a4)

    def _update_property_keys(self, new_frame=None):
        """
        Update the internal list of molecular properties that can be sliced. If a frame is passed,
        update from that frame's properties. Otherwise, update from the entire list of stored frames
        """
        if self._property_keys is None:
            self._property_keys = set()

        if new_frame is None:
            for frame in self:
                self._property_keys.update(frame.keys())
        else:
            self._property_keys.update(new_frame)

    def __dir__(self):
        return list(self._property_keys.union(dir(self.__class__)).union(self.__dict__))

    def __getitem__(self, item):
        return self.frames[item]

    def __getattr__(self, item):
        # TODO: move slicing to frames, so that this will work with __getitem__ as well
        # TODO: prevent identical recursion (so __getattr__(foo) can't call __getattr__(foo))
        if self._property_keys is None:
            self._update_property_keys()

        if item in self._property_keys:  # return a list of properties for each frame
            return self.slice_frames(item)
        else:
            raise AttributeError('Frame %s has no attribute %s' % (self, item))

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

    def slice_frames(self, key, missing=None):
        """ Return an array of giving the value of ``key`` at each frame.

        Args:
            key (str): name of the property, e.g., time, potential_energy, annotation, etc
            missing: value to return if a given frame does not have this property

        Returns:
            moldesign.units.Vector: vector containing the value at each frame, or the value given
                in the ``missing`` keyword) (len= `len(self)` )
        """
        has_units = True
        result = []
        for f in self.frames:
            val = f.get(key, None)
            if not issubclass(type(val), u.MdtQuantity):
                has_units = False
            result.append(val)
        if has_units:
            result = u.array([frame.get(key, None) for frame in self.frames])
            return u.default.convert(result)
        else:
            return np.array(result)

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

    DONOTAPPLY = set(['kinetic_energy'])

    def apply_frame(self, frame):
        """
        Reconstruct the underlying molecule with the given frame.
        Right now, any data not passed is ignored, which may result in properties that aren't synced up
        with each other ...
        """
        # TODO: need to prevent multiple things using the _tempmol from conflicting with each other
        m = self._tempmol
        for prop in self.MOL_ATTRIBUTES:
            if prop in self.DONOTAPPLY:
                continue
            if prop in frame:
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

    def plot(self, x, y, **kwargs):
        """ Create a matplotlib plot of property x against property y

        Args:
            x,y (str): names of the properties
            **kwargs (dict): kwargs for :meth:`matplotlib.pylab.plot`

        Returns:
            List[matplotlib.lines.Lines2D]: the lines that were plotted

        """
        from matplotlib import pylab
        xl = yl = None
        if type(x) is str:
            strx = x
            x = getattr(self, x)
            xl = '%s / %s' % (strx, x.units)
        if type(y) is str:
            stry = y
            y = getattr(self, y)
            yl = '%s / %s' % (stry, y.units)
        plt = pylab.plot(x, y, **kwargs)
        pylab.xlabel(xl); pylab.ylabel(yl); pylab.grid()
        return plt

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
        elif type(reference_frame) == int:
            iframe = reference_frame
            relative_alignment = True
        else:
            relative_alignment = False

        if relative_alignment:
            for i in xrange(iframe+1, self.num_frames):
                self.frames[i].wfn.align_orbital_phases(
                    self.frames[i-1].wfn)
            for i in xrange(iframe-1, -1, -1):
                self.frames[i].wfn.align_orbital_phases(
                    self.frames[i+1].wfn)
        else:
            for frame in self.frames:
                frame.wfn.align_orbital_phases(reference_frame.wfn)
