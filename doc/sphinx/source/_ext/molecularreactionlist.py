# AUTO-GENERATE DOCUMENTATION FOR DREAM MOLECULAR REACTIONS
#
# Builds a table of the reactions defined in
# py/DREAM/Settings/MolecularReactionData.py, and of the species in
# SPECIES_Z, so that the documentation
# always matches the reactions available in the Python interface.
#
# Usage:
#
#   .. molecularspecieslist::
#
#   .. molecularreactionlist::
#      :classification: CHARGE_EXCHANGE
#
# If ':classification:' is omitted, all reactions are listed.
# ##########################

import re

from docutils import nodes
from docutils.parsers.rst import Directive, directives
from docutils.parsers.rst.directives.tables import ListTable

from DREAM.Settings import MolecularReactionData


class MolecularReactionList(ListTable, Directive):

    option_spec = {
        'classification': directives.unchanged
    }


    def run(self):
        """
        Generate the list of molecular reactions.
        """
        classification = self.options.get('classification', None)

        reactions = MolecularReactionData.MOLECULAR_REACTIONS
        if classification is not None:
            reactions = [r for r in reactions if r['classification'] == classification]

        col_widths = [3, 4, 2]
        header_rows = 1
        stub_columns = 0

        table_data = []
        # Add header
        table_data.append([
            nodes.paragraph(text='Name'),
            nodes.paragraph(text='Reaction'),
            nodes.paragraph(text='Classification')
        ])

        # Add rows to table
        for r in reactions:
            reaction = self.format_side(r['reactants']) + r' \rightarrow ' + self.format_side(r['products'])
            table_data.append([
                nodes.paragraph('', '', nodes.literal(text=r['name'])),
                nodes.paragraph('', '', nodes.math(reaction, reaction)),
                nodes.paragraph('', '', nodes.literal(text=r['classification']))
            ])

        # Construct the table node
        table_node = self.build_table_from_list(
            table_data, col_widths, header_rows, stub_columns
        )

        return [table_node]


    def format_side(self, participants):
        """
        Convert a list of (name, Z0, coefficient) tuples into LaTeX,
        e.g. [('e', -1, 1), ('D2', 0, 1)] -> 'e^- + D_2'.
        """
        terms = []
        for name, Z0, coefficient in participants:
            s = self.format_species(name, Z0)
            if coefficient > 1:
                s = '{}\\,{}'.format(coefficient, s)
            terms.append(s)

        return ' + '.join(terms)


    def format_species(self, name, Z0):
        """
        Convert a species name and charge state into LaTeX,
        e.g. ('D2', 1) -> '\\mathrm{D}_2^{+}'.
        """
        if name == 'e':
            return r'\mathrm{e}^-'
        elif name == 'RE':
            return r'\mathrm{e}^-_\mathrm{RE}'

        # Write the number of atoms as subscripts (D2 -> D_2)
        s = re.sub(r'([A-Za-z]+)(\d+)', r'\\mathrm{\1}_{\2}', name)
        if s == name:
            s = '\\mathrm{{{}}}'.format(name)

        if Z0 == 1:
            s += '^{+}'
        elif Z0 > 1:
            s += '^{{{}+}}'.format(Z0)

        return s


class MolecularSpeciesList(ListTable, Directive):

    def run(self):
        """
        Generate the list of species appearing in the molecular reactions.
        """
        col_widths = [1, 1, 1]
        header_rows = 1
        stub_columns = 0

        table_data = []
        # Add header
        table_data.append([
            nodes.paragraph(text='Name'),
            nodes.paragraph(text='Z'),
            nodes.paragraph(text='Type')
        ])

        # Add rows to table
        for name, Z in MolecularReactionData.SPECIES_Z.items():
            # Count the atoms in the name (D2 -> 2, ArD -> 2, He -> 1)
            nAtoms = sum(int(n) if n else 1 for _, n in re.findall(r'([A-Z][a-z]?)(\d*)', name))
            table_data.append([
                nodes.paragraph('', '', nodes.literal(text=name)),
                nodes.paragraph(text=str(Z)),
                nodes.paragraph(text='Molecule' if nAtoms > 1 else 'Atom')
            ])

        # Construct the table node
        table_node = self.build_table_from_list(
            table_data, col_widths, header_rows, stub_columns
        )

        return [table_node]


def setup(app):
    app.add_directive('molecularreactionlist', MolecularReactionList)
    app.add_directive('molecularspecieslist', MolecularSpeciesList)

    return {
        'version': '0.1',
        'parallel_read_safe': True,
        'parallel_write_safe': True
    }
