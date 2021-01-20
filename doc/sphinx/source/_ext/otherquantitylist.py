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
        'source': directives.path,
        'groups': directives.unchanged
    }


    def run(self):
        """
        Generate the list of other quantities.
        """
        source = self.options.get('source', '../../../../../src/OtherQuantityHandler.cpp')
        groups = self.options.get('groups', 'no').lower() == 'yes'

        table_data = []

        #col_widths = self.get_column_widths(2)
        col_widths = [1, 2]
        header_rows = 1
        stub_columns = 0

        table_data = []
        # Add header
        if groups:
            table_data.append([nodes.paragraph(text='Name'), nodes.paragraph(text='Included quantities')])
        else:
            table_data.append([nodes.paragraph(text='Name'), nodes.paragraph(text='Description')])

        if groups:
            qd = self.load_groups(source)
        else:
            qd = self.load_quantities(source)

        # Add rows to table
        for row in qd:
            if type(row[1]) == str:
                table_data.append([nodes.literal(text=row[0]), nodes.paragraph(text=row[1])])
            else:
                table_data.append([nodes.literal(text=row[0]), row[1]])

        # Construct the table node
        table_node = self.build_table_from_list(
            table_data, col_widths, header_rows, stub_columns
        )

        return [table_node]


    def load_string(self, s):
        """
        Load the next string (i.e. text in quotation marks) from
        the given line 's'.
        """
        n = len(s)

        # Locate beginning quotation mark
        i = 0
        while i < n and s[i] != '"': i += 1

        if i == n:
            return None, None

        # Load string
        i += 1
        j = 0

        while i+j < n and s[i+j] != '"' and s[i+j-1] != '\\':
            j += 1

        S = s[i:(i+j)]

        return S, i+j+1


    def load_groups(self, source):
        """
        Load the list of other quantity groups and their descriptions from the
        DREAM source code.
        """
        with open(os.path.abspath(source), 'r') as f:
            src = f.readlines()

        grps = []
        match = 'this->groups['
        I = 0
        nLines = len(src)
        while I < nLines:
            l = src[I].strip()

            if l.startswith(match):
                i = len(match)
                group, _ = self.load_string(l[i:])

                # Load quantities included in the group
                q = ''
                if l[-1] == '{':
                    q = []
                    I += 1
                    while I < nLines:
                        l = src[I].strip()

                        # End of list of quantities?
                        if l.startswith('};'):
                            break

                        # Otherwise, load the list of quantities
                        i = 0
                        s, j = self.load_string(l)
                        while s is not None:
                            #ln = nodes.line()
                            #ln.append(nodes.literal(text=s))
                            q.append(s)
                            i += j

                            s, j = self.load_string(l[i:])

                        I += 1

                    # Sort q
                    q.sort(key=lambda x : x.lower())

                    # Convert q into list of lines
                    Q = []
                    for quant in q:
                        ln = nodes.line()
                        ln.append(nodes.literal(text=quant))
                        Q.append(ln)

                    #grps.append([group, ', '.join(q)])
                    #grps.append([group, pg])
                    grps.append([group, Q])
                else:
                    q = "All quantities starting with ``{}/``.".format(group)

            I += 1

        return grps


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
                i = len('DEF_')
                name, j = self.load_string(l[i:])
                desc, _ = self.load_string(l[(i+j):])
                """
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
                """

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

