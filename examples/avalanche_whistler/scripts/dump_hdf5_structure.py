#!/usr/bin/env python3
#
# Dump the structure of a DREAM HDF5 output file as an indented text tree.
#
# Usage:
#   python dump_hdf5_structure.py <input.h5>                     # print to stdout
#   python dump_hdf5_structure.py <input.h5> -o structure.txt    # write to file
#

import argparse
import os

import h5py


def dump_group(group, lines, depth):
    """Recursively append '- name (Group|Dataset)' lines for every member."""
    for name, obj in group.items():
        kind = 'Group' if isinstance(obj, h5py.Group) else 'Dataset'
        lines.append('  ' * depth + f'- {name} ({kind})')
        if isinstance(obj, h5py.Group):
            dump_group(obj, lines, depth + 1)


def main():
    parser = argparse.ArgumentParser(
        description='Dump the structure of an HDF5 file as an indented text tree')
    parser.add_argument('input', help='HDF5 file to inspect')
    parser.add_argument('-o', '--output', type=str, default=None,
                        help='Output text file (default: print to stdout)')
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"Error: file not found: {args.input}")
        return

    with h5py.File(args.input, 'r') as f:
        lines = []
        dump_group(f, lines, depth=0)

    text = '\n'.join([
        '',
        f'HDF5 File Structure: {os.path.basename(args.input)}',
        '=' * 50,
    ] + lines) + '\n'

    if args.output:
        with open(args.output, 'w') as fh:
            fh.write(text)
        print(f"Structure written to {args.output}")
    else:
        print(text, end='')


if __name__ == '__main__':
    main()
