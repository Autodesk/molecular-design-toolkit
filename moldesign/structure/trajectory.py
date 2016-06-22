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

import ipywidgets as ipy
import numpy as np

import moldesign as mdt
from moldesign import units as u
from moldesign.units import default
from moldesign.utils import DotDict
from moldesign.widgets import ui

__all__ = 'Trajectory TrajectoryViz'.split()


class SubSelection(object):
    """
    Descriptors to get bits of the trajectory
    trajectory.atoms[3].position -> array of positions
    trajectory.atoms[5].distance( mol.chains['B'].residues[5] ) -> array of distances
    trajectory.chains['A'].residue[23].com -> array of COMs
    NOT IMPLEMENTED YET
    """


class Frame(DotDict):
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


class Trajectory(object):
    """
    A ``Trajectory`` stores information about a molecule's motion and how its properties change as
    it moves. A trajectory object contains 1) a reference to the :class:`Molecule` it describes, and
    2) a list of :class:`Frame` objects, each one containing a snapshot of the molecule at a
    particular point in its motion.
    """

    # TODO: the current implementation does not support cases where the molecular topology changes
    # TODO: allow caching to disk for very large trajectories

    MOL_ATTRIBUTES = ['positions', 'momenta', 'time']
    """List[str]: Always store these molecular attributes"""

    def __getstate__(self):
        return self.__dict__.copy()

    def __setstate__(self, d):
        self.__dict__.update(d)

    def __init__(self, mol, unit_system=default):
        """ Initialization: create a new, empty trajectory based on the passed molecule

        Args:
            mol (moldesign.Molecule): the trajectory will describe the motion of this molecule
            unit_system (u.UnitSystem): convert all attributes to this unit system
        """
        self._init = True
        self.info = "Trajectory"
        self.frames = []
        self.mol = mol
        self.unit_system = unit_system
        self._property_keys = None
        self._tempmol = mdt.Molecule(self.mol.atoms, copy_atoms=True)
        self._tempmol.dynamic_dof = self.mol.dynamic_dof
        self._viz = None
        self.atoms = [TrajObject(a, self) for a in mol.atoms]
        self.chains = [TrajObject(c, self) for c in mol.chains]
        self.residues = [TrajObject(r, self) for r in mol.residues]

    """ Attributes:
        mol (moldesign.Molecule): the molecule object that this trajectory comes from
        frames (List[Frame]): a list of the trajectory frames in the order they were created
        info (str): text describing this trajectory
        unit_system (u.UnitSystem): convert all attributes to this unit system
    """

    @property
    def num_frames(self):
        """int: number of frames in this trajectory"""
        return len(self)

    def __len__(self):
        """overrides len(trajectory) to return number of frames"""
        return len(self.frames)

    def draw3d(self):
        """TrajectoryViz: create a trajectory visualization"""
        self._viz = TrajectoryViz(self)
        return self._viz
    draw = draw3d  # synonym for backwards compatibility

    def draw_orbitals(self, align=True):
        """TrajectoryOrbViz: create a trajectory visualization"""
        if align: self.align_orbital_phases()
        self._viz = TrajectoryOrbViz(self)
        return self._viz

    def __str__(self):
        return 'Trajectory for molecule "%s" (%d frames)' % (self.mol, self.num_frames)

    def __repr__(self):
        try:
            return '<%s>' % str(self)
        except Exception:
            return '<Trajectory object @ %s (exception in repr)>' % hex(id(self))

    def __add__(self, other):
        newtraj = Trajectory(self, unit_system=self.unit_system)
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
        """
        Two possibilities here:
        If "item" is one of the stored properties, return an frame slice
        If "item" is a molecular entity, return a trajectory with positions, momenta, and forces that
        just refer to that entity.
        """
        # TODO: prevent identical recursion (so __getattr__(foo) can't call __getattr__(foo))
        if self._property_keys is None:
            self._update_property_keys()

        if item in self._property_keys:  # return a list of properties for each frame
            return self.slice_frames(item)
        else:
            raise AttributeError

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
        indices = []
        for atom in atoms:
            indices.extend(range(*atom.parent_slice.indices(self.mol.ndims)))
        indices = np.array(indices)

        rmsds = []
        for f in self.frames:
            diff = (refpos[indices] - f.positions[indices])
            rmsds.append(np.sqrt(diff.dot(diff)/len(atoms)))
        return u.to_units_array(rmsds).defunits()

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
            if not issubclass(type(val), u.BuckyballQuantity):
                has_units = False
            result.append(val)
        if has_units:
            result = u.to_units_array([frame.get(key, None) for frame in self.frames])
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
                    moldesign.core.helpers.kinetic_energy(frame.momenta, self.mol.dim_masses))
            else:
                convert_units = False
                energies.append(None)
        if convert_units:
            arr = u.to_units_array(energies)
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
                temps.append(moldesign.core.helpers.kinetic_temperature(energy, dof))
            else:
                convert_units = False
                temps.append(None)
        if convert_units:
            arr = u.to_units_array(temps)
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
        m.properties = moldesign.structure.molecule.MolecularProperties(m)
        for attr in frame:
            m.properties[attr] = frame[attr]

    def write_trajectory(self, filename=None, format=None, overwrite=True):
        """
        Write to a file (if filename provided) or file-like buffer

        Args:
            filename (str): name of file (return a file-like object if not passed)
            format (str): file format (guessed from filename if None)
            overwrite (bool): overwrite filename if it exists

        Returns:
            StringIO: file-like object (only if filename not passed)
        """
        tempmol = self._tempmol
        if filename and (not overwrite) and os.path.exists(filename):
            raise IOError('%s exists' % filename)
        if not filename:
            fileobj = StringIO()
        else:
            fileobj = open(filename, 'w')

        if format is None:
            format = filename.split('.')[-1]

        for frame in self.frames:
            self.apply_frame(frame)
            fileobj.write(tempmol.write(format=format))

        if filename is None:
            fileobj.seek(0)
            return fileobj
        else:
            fileobj.close()

    def plot(self, x, y, **kwargs):
        """ Create a matplotlib plot of property x against property y

        Args:
            x,y (str): names of the properties
            **kwargs (dict): kwargs for :function:`matplotlib.pylab.plot`

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
                self.frames[i].electronic_state.align_orbital_phases(
                    self.frames[i-1].electronic_state)
            for i in xrange(iframe-1, -1, -1):
                self.frames[i].electronic_state.align_orbital_phases(
                    self.frames[i+1].electronic_state)
        else:
            for frame in self.frames:
                frame.electronic_state.align_orbital_phases(reference_frame.electronic_state)


class TrajectoryViz(ui.SelectionGroup):
    def __init__(self, trajectory, **kwargs):
        """
        :type trajectory: moldesign.trajectory.Trajectory
        :param kwargs: kwargs for widget
        :return:
        """
        self.default_fps = 10
        self.traj = trajectory
        self.pane = ipy.VBox()
        trajectory.apply_frame(trajectory.frames[0])
        self.viewer, self.view_container = self.make_viewer()
        for frame in self.traj.frames[1:]:
            self.viewer.append_frame(positions=frame.positions.reshape((trajectory.mol.num_atoms, 3)),
                                     render=False)
        self.make_controls()
        self.pane.children = [self.view_container, self.controls]
        super(TrajectoryViz, self).__init__([self.pane, ui.AtomInspector()], **kwargs)
        self.update_selections('initialization', {'framenum': 0})

    def make_viewer(self):
        viewer = self.traj._tempmol.draw3d(style='licorice')
        viewer.show_unbonded()
        return viewer, viewer

    def make_controls(self):
        controls = []
        self.play_button = ipy.Button(description=u"\u25B6")
        self.text_view = FrameInspector(self.traj)
        controls.append(self.text_view)
        self.play_button.on_click(self._play)
        self.slider = ui.create_value_selector(ipy.IntSlider, value_selects='framenum',
                                               value=0, description='Frame:',
                                               min=0, max=len(self.traj)-1)
        self.playbox = ipy.HBox([self.play_button, self.slider])
        controls.append(self.playbox)
        self.controls = ipy.VBox(controls)

    def _play(self, *args, **kwargs):
        self.animate()

    def animate(self, fps=None):
        fps = mdt.utils.if_not_none(fps, self.default_fps)
        self.slider.value = 0
        spf = 1.0 / fps
        for i in xrange(self.viewer.num_frames):
            t0 = time.time()
            self.slider.value = i
            dt = time.time() - t0
            if dt < spf:
                time.sleep(spf-dt)

    def __getattr__(self, item):
        """Users can run viz commands directly on the trajectory,
        e.g., trajectory.cpk(atoms=mol.chains['B'])"""
        # TODO: modify __dir__ to match
        return getattr(self.viewer, item)


class TrajectoryOrbViz(TrajectoryViz):
    def make_viewer(self):
        viewframe = self.traj._tempmol.draw_orbitals()
        viewframe.viewer.frame_change_callback = self.on_frame_change
        return viewframe.viewer, viewframe

    def on_frame_change(self, framenum):
        self.traj.apply_frame(self.traj.frames[framenum])
        self.viewer.wfn = self.traj._tempmol.electronic_state


class FrameInspector(ipy.HTML, ui.Selector):
    def __init__(self, traj, **kwargs):
        self.traj = traj
        super(FrameInspector, self).__init__(**kwargs)

    def handle_selection_event(self, selection):
        if 'framenum' not in selection:
            return
        else:
            framenum = selection['framenum']

        if hasattr(self.traj[framenum], 'time') and self.traj[framenum].time is not None:
            result = 'Time: %s; ' % self.traj[framenum].time.defunits()
        else:
            result = ''

        try:
            result += self.traj[framenum].annotation
        except (KeyError, AttributeError):
            pass

        self.value = result


class TrajObject(object):  # under construction ...
    def __init__(self, obj, traj):
        """
        :param obj:
        :type traj: Trajectory
        :return:
        """
        self.obj = obj
        self.traj = traj

    def distance(self, other):
        retvals = []
        for frame in self.traj.frames:
            self.traj.apply_frame(frame)
            retvals.append(self.obj.distance(other.obj))
        return retvals
