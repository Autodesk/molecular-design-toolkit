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

from moldesign import utils
from moldesign.helpers import colormap


class ColorMixin(object):
    def color_by(self, atom_callback, atoms=None, mplmap='auto', render=True,
                 force_cmap=False):
        """
        Color atoms according to either:
          * an atomic attribute (e.g., 'chain', 'residue', 'mass')
          * a callback function that accepts an atom and returns a color or a category

        Args:
            atom_callback (callable OR str): callable f(atom) returns color OR
                category OR an atom attribute (e.g., ``atnum, mass, residue.type``)
            atoms (moldesign.molecules.AtomContainer): atoms to color (default: self.atoms)
            mplmap (str): name of the matplotlib colormap to use if colors aren't explicitly
               specified)
            force_cmap (bool): force the use of a colormap
            render (bool): draw these changes immediately

        Notes:
            If you'd like to explicitly specify colors, the callback can return color
            specifications as an HTML string (``'#1234AB'``), a hexadecimal integer (
            ``0x12345AB``), or a CSS3 color keyword (``'green'``, ``'purple'``, etc., see
            https://developer.mozilla.org/en-US/docs/Web/CSS/color_value)

            If the callback returns an integer, it may be interpreted as a color spec (since RGB
            colors are just hexadecimal integers). Use ``force_cmap=True`` to force the creation
            of a colormap.

        Returns:
            dict: mapping of categories to colors
        """
        atoms = utils.if_not_none(atoms, self.mol.atoms)
        if isinstance(atom_callback, basestring):
            # shortcut to use strings to access atom attributes, i.e. "ff.partial_charge"
            attrs = atom_callback.split('.')

            # make sure that whatever value is returned doesn't get interpreted as a color
            force_cmap = True

            def atom_callback(atom):
                obj = atom
                for attr in attrs:
                    obj = getattr(obj, attr)
                return obj

        colors = utils.Categorizer(atom_callback, atoms)

        if force_cmap:
            name_is_color = [False]
        else:
            name_is_color = map(utils.is_color, colors.keys())

        if len(colors) <= 1:
            colors = {'gray': atoms}

        elif not all(name_is_color):
            assert not any(name_is_color), \
                "callback function returned a mix of colors and categories"
            categories = colors
            cats = categories.keys()
            # If there are >256 categories, this is a many-to-one mapping
            colornames = colormap(cats, mplmap=mplmap)
            colors = {c: [] for c in colornames}
            for cat, color in zip(cats, colornames):
                colors[color].extend(categories[cat])

        self.set_colors(colors, render=render)

