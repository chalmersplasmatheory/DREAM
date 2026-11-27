# AUTO-GENERATE DOCUMENTATION FOR DREAM STATIC QUANTITIES
#
# Written by: Mathias Hoppe
# ##########################

from docutils import nodes
from docutils.parsers.rst import Directive, directives
from docutils.parsers.rst.directives.tables import ListTable

import os


class StaticQuantityList(ListTable, Directive):

    option_spec = {
        'source': directives.path
    }


    def run(self):
        """
        Generate the list of static quantities.
        """
        source = self.options.get('source', '../../../../../src/StaticQuantityHandler.cpp')

        col_widths = [1, 2]
        header_rows = 1
        stub_columns = 0

        # Add header
        table_data = [[nodes.paragraph(text='Name'), nodes.paragraph(text='Description')]]

        # Add rows to table
        for name, desc in self.load_quantities(source):
            table_data.append([nodes.literal(text=name), nodes.paragraph(text=desc)])

        # Construct the table node
        table_node = self.build_table_from_list(
            table_data, col_widths, header_rows, stub_columns
        )

        return [table_node]


    def load_string(self, s):
        """
        Load the next string (i.e. text in quotation marks) from
        the given line 's'. Returns the string and the index of the
        character immediately following the closing quotation mark.
        """
        n = len(s)

        # Locate beginning quotation mark
        i = s.find('"')
        if i < 0:
            return None, None

        # Load string (respecting escaped quotation marks)
        i += 1
        j = i
        while j < n and not (s[j] == '"' and s[j-1] != '\\'):
            j += 1

        S = s[i:j].replace('\\"', '"')

        return S, j+1


    def load_quantities(self, source):
        """
        Load the list of static quantities and descriptions from the DREAM
        source code.
        """
        with open(os.path.abspath(source), 'r') as f:
            src = f.readlines()

        quants = []
        for line in src:
            l = line.strip()

            # Does this line contain a StaticQuantity definition?
            # (the '#define DEF_QTY' line is skipped automatically)
            if l.startswith('DEF_QTY('):
                i = len('DEF_QTY(')
                name, j = self.load_string(l[i:])
                if name is None:
                    continue

                desc, _ = self.load_string(l[(i+j):])
                if desc is None:
                    desc = ''

                quants.append([name, desc])

        # Sort list (based on name)
        quants.sort(key=lambda x : x[0].lower())
        return quants


def setup(app):
    app.add_directive("staticquantitylist", StaticQuantityList)

    return {
        'version': '1.0',
        'parallel_read_safe': True,
        'parallel_write_safe': True
    }
