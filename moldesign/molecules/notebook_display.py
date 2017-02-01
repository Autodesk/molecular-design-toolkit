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


def sequence_view(chain):
    """ Creates a simple primary structure viewer
    """
    from IPython.display import HTML

    html = ['<span style="word-wrap:break-word">']
    for residue in chain.residues:
        html.append('<span title="%s (idx %d), chain %s">%s</a>' %
                    (residue.name, residue.index, residue.chain.name,
                     residue.code))
    html.append("</span>")

    return HTML(''.join(html))

class MolNotebookMixin(object):
    """ Methods for displaying molecular information in Jupyter notebooks

    Note:
        This is a mixin class designed only to be mixed into the :class:`Molecule` class. Routines
        are separated are here for code organization only - they could be included in the main
        Molecule class without changing any functionality
    """

    def markdown_summary(self):
        """A markdown description of this molecule.

        Returns:
            str: Markdown"""
        # TODO: remove leading underscores for descriptor-protected attributes
        lines = ['### Molecule: "%s" (%d atoms)' % (self.name, self.natoms)]

        description = self.metadata.get('description', None)
        if description is not None:
            description = self.metadata.description[:5000]
            url = self.metadata.get('url', None)
            if url is not None:
                description = '<a href="%s" target="_blank">%s</a>' % \
                              (url, description)
            lines.append(description)

        lines.extend([
                 '**Mass**: {:.2f}'.format(self.mass),
                 '**Formula**: %s' % self.get_stoichiometry(html=True),
                 '**Charge**: %s' % self.charge])

        if self.energy_model:
            lines.append('**Potential model**: %s' % str(self.energy_model))

        if self.integrator:
            lines.append('**Integrator**: %s' % str(self.integrator))

        if self.is_biomolecule:
            lines.extend(self.biomol_summary_markdown())

        return '\n\n'.join(lines)

    def _repr_markdown_(self):
        return self.markdown_summary()

    def biomol_summary_markdown(self):
        """A markdown description of biomolecular structure.

        Returns:
            str: Markdown string"""
        lines = []
        if len(self.residues) > 1:
            table = self.get_residue_table()
            lines.append('### Residues')
            # extra '|' here may be workaround for a bug in ipy.markdown?
            lines.append(table.markdown(replace={0: ' '}) + '|')

            lines.append('### Biopolymer chains')
            seqs = []
            for chain in self.chains:
                seq = chain._get_sequence(_html=True)
                if not seq.strip():  # don't write anything if there's no sequence
                    continue

                seqs.append('**%s**: %s' % (chain.name, seq))
            lines.append('<br>'.join(seqs))

        return lines

    def get_residue_table(self):
        """Creates a data table summarizing this molecule's primary structure.

        Returns:
            moldesign.utils.MarkdownTable"""
        table = utils.MarkdownTable(*(['chain'] +
                                      'protein dna rna unknown water solvent'.split()))
        for chain in self.chains:
            counts = {}
            unk = []
            for residue in chain.residues:
                cat = residue.type
                if cat == 'unknown':
                    unk.append(residue.name)
                counts[cat] = counts.get(cat, 0) + 1
            counts['chain'] = '<b>%s</b>' % chain.name
            if 0 < len(unk) <= 4:
                counts['unknown'] = ','.join(unk)
            table.add_line(counts)
        return table

    def get_stoichiometry(self, html=False):
        """ Return this molecule's stoichiometry

        Returns:
            str
        """
        counts = {}
        for atom in self.atoms:
            counts[atom.symbol] = counts.get(atom.symbol, 0) + 1

        my_elements = sorted(counts.keys())
        if html: template = '%s<sub>%d</sub>'
        else: template = '%s%d'
        return ''.join([template % (k, counts[k]) for k in my_elements])