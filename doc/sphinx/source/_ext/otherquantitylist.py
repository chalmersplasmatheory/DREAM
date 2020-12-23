# AUTO-GENERATE DOCUMENTATION FOR DREAM OTHER QUANTITIES
#
# Written by: Mathias Hoppe
# ##########################

from docutils import nodes
from docutils.parsers.rst import Directive, directives
from docutils.parsers.rst.directives.tables import ListTable

import os

from sphinx.locale import _
from sphinx.util.docutils import SphinxDirective


class OtherQuantityList(ListTable, Directive):
    
    option_spec = {
        'source': directives.path
    }


    def run(self):
        """
        Generate the list of other quantities.
        """
        source = self.options.get('source', '../../../../../src/OtherQuantityHandler.cpp')

        table_data = []

        #col_widths = self.get_column_widths(2)
        col_widths = [1, 2]
        header_rows = 1
        stub_columns = 0

        table_data = []
        # Add header
        table_data.append([nodes.paragraph(text='Name'), nodes.paragraph(text='Description')])

        qd = self.load_quantities(source)

        # Add rows to table
        for row in qd:
            table_data.append([nodes.literal(text=row[0]), nodes.paragraph(text=row[1])])

        # Construct the table node
        table_node = self.build_table_from_list(
            table_data, col_widths, header_rows, stub_columns
        )

        return [table_node]


    def load_quantities(self, source):
        """
        Load the list of other quantities and descriptions from the DREAM
        source code.
        """
        with open(os.path.abspath(source), 'r') as f:
            src = f.readlines()

        quants = []
        for line in src:
            l = line.strip()
            
            # Does this line contain an OtherQuantity definition?
            if l.startswith('DEF_'):
                i = 4
                while l[i] != '"': i += 1

                # Load quantity name
                i += 1
                j = 0
                while l[i+j] != '"' and l[i+j-1] != '\\':
                    j += 1

                name = l[i:(i+j)]

                # Locate quantity description
                i += j+1
                while l[i] != '"':
                    i += 1

                # Load quantity description
                i += 1
                j = 0
                while l[i+j] != '"' and l[i+j-1] != '\\':
                    j += 1

                desc = l[i:(i+j)]

                quants.append([name, desc])

        # Sort list (based on name)
        quants.sort(key=lambda x : x[0].lower())
        return quants


def setup(app):
    app.add_directive("otherquantitylist", OtherQuantityList)

    return {
        'version': '1.0',
        'parallel_read_safe': True,
        'parallel_write_safe': True
    }

