# Simple wrapper for running 'dreami'

import os
import pathlib
import subprocess
import tempfile

from . DREAMException import DREAMException
from . DREAMOutput import DREAMOutput
from . DREAMSettings import DREAMSettings


DREAMPATH = None

def locatedream():
    global DREAMPATH

    try:
        DREAMPATH = os.environ['DREAMPATH']

        if DREAMPATH[-1] == '/':
            DREAMPATH = DREAMPATH[:-1]

    except KeyError: pass

    if DREAMPATH is None:
        DREAMPATH = (pathlib.Path(__file__).parent / '..' / '..').resolve().absolute()

        if not os.path.isfile('{}/build/iface/dreami'.format(DREAMPATH)):
            raise DREAMException("Unable to locate the DREAMi executable. Try to set the 'DREAMPATH' environment variable.")


def runiface(settings, outfile=None):
    """
    Run 'dreami' with the specified settings (which may be either
    a 'DREAMSettings' object or the name of a file containing the
    settings).

    settings: 'DREAMSettings' object or name of file containing settings.
    outfile:  Name of file to write output to (default: 'output.h5')
    """
    global DREAMPATH

    infile = None
    if isinstance(settings, DREAMSettings):
        infile = next(tempfile._get_candidate_names())+'.h5'
        settings.save(infile)
    else:
        infile = settings

    deleteOutput = False
    if outfile is None:
        deleteOutput = True
        outfile = next(tempfile._get_candidate_names())+'.h5'

    p = subprocess.Popen(['{}/build/iface/dreami'.format(DREAMPATH), infile, '-o', outfile], stderr=subprocess.PIPE)
    stderr_data = p.communicate()[1].decode('utf-8')

    os.remove(infile)

    if p.returncode != 0:
        print(stderr_data)
        raise DREAMException("DREAMi exited with a non-zero exit code: {}".format(p.returncode))

    obj = DREAMOutput(outfile)

    if deleteOutput:
        os.remove(outfile)

    return obj

locatedream()

