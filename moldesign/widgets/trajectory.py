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
import time

import ipywidgets as ipy

import moldesign as mdt

from moldesign.uibase.components import AtomInspector
from moldesign.uibase import selector


class TrajectoryViewer(selector.SelectionGroup):
    """ 3D representation, with animation controls, for a trajectory.

    Users will typically instantiate this using ``trajectory.draw()``

    Args:
        display (bool): immediately display this to the notebook (default: False)
        **kwargs (dict): keyword arguments for :class:`ipywidgets.Box`
    """

    def __init__(self, trajectory, display=False, **kwargs):
        from IPython.display import display as displaynow

        self.default_fps = 10
        self.traj = trajectory
        self.pane = ipy.VBox()
        trajectory.apply_frame(trajectory.frames[0])
        self.viewer, self.view_container = self.make_viewer()
        for frame in self.traj.frames[1:]:
            self.viewer.append_frame(positions=frame.positions,
                                     wfn=frame.get('wfn', None),
                                     render=False)
        self.make_controls()
        self.pane.children = [self.view_container, self.controls]
        super(TrajectoryViewer, self).__init__([self.pane, AtomInspector()], **kwargs)
        self.update_selections('initialization', {'framenum': 0})
        if display:
            displaynow(self)

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
        self.slider = selector.create_value_selector(ipy.IntSlider,
                                                     value_selects='framenum',
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


class TrajectoryOrbViewer(TrajectoryViewer):
    def make_viewer(self):
        viewframe = self.traj._tempmol.draw_orbitals()
        return viewframe.viewer, viewframe


class FrameInspector(ipy.HTML, selector.Selector):
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
