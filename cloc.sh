#!/usr/bin/env bash
#
# This script runs the 'cloc' command (for "Count Lines Of Code")
# in a manner appropriate for counting the number of lines of
# code in DREAM.
#

# Ensure that 'cloc' is installed first
if ! which cloc &>/dev/null; then
    echo "ERROR: CLOC is not installed. Please install it first."
    echo "(e.g. 'sudo apt install cloc')"
    exit 1
fi

directories="fvm iface include src"
if [ $# -eq 1 ]; then
    if [ "$1" == "--conservative" ]; then
        directories="$directories"
    elif [ "$1" == "--liberal" ]; then
        directories="$directories py tools tests"
    elif [ "$1" == "--help" ]; then
        echo -e "./cloc.sh [--conservative] [--liberal]"
        echo -e "   COUNT LINES OF CODE IN DREAM\n"

        echo "FLAGS"
        echo -e "  --conservative      Estimate number of lines conservatively (counting"
        echo -e "                      only code \"that counts\")"
        echo -e "  --help              Show this help message"
        echo -e "  --liberal           Estimate number of lines liberally (counting most"
        echo -e "                      everything in the code base)"
    else
        echo "ERROR: Unrecognized flag '$1'. Use '--help' for usage info."
        exit 1
    fi
else
    directories="$directories py"
fi

cloc $directories --fullpath --not-match-d="src/Atomics"

