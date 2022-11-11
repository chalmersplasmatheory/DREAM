# Simple wrapper for running 'dreami'

import os
import pathlib
import subprocess
import tempfile

from . DREAMException import DREAMException
from . DREAMOutput import DREAMOutput
from . DREAMSettings import DREAMSettings
from subprocess import TimeoutExpired

from subprocess import TimeoutExpired

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
            #raise DREAMException("Unable to locate the DREAMi executable. Try to set the 'DREAMPATH' environment variable.")
            print("WARNING: Unable to locate the DREAMi executable. Try to set the 'DREAMPATH' environment variable.")


def runiface(settings, outfile=None, quiet=False, timeout=None):
    """
    Run 'dreami' with the specified settings (which may be either
    a 'DREAMSettings' object or the name of a file containing the
    settings).

    settings: 'DREAMSettings' object or name of file containing settings.
    outfile:  Name of file to write output to (default: 'output.h5')
    """
    global DREAMPATH

    deleteOutput = False
    if outfile is None:
        deleteOutput = True
        outfile = next(tempfile._get_candidate_names())+'.h5'

    infile = None
    if isinstance(settings, DREAMSettings):
        infile = next(tempfile._get_candidate_names())+'.h5'
        settings.output.setFilename(outfile)
        settings.save(infile)
    else:
        infile = settings

    errorOnExit = 0
    p = None
    obj = None
    stderr_data = None
    try:
        if quiet:
            p = subprocess.Popen(['{}/build/iface/dreami'.format(DREAMPATH), infile], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        else:
            p = subprocess.Popen(['{}/build/iface/dreami'.format(DREAMPATH), infile], stderr=subprocess.PIPE)

        try:
            stderr_data = p.communicate(timeout=timeout)[1].decode('utf-8')

            if p.returncode != 0:
                errorOnExit = 1
            else:
                obj = DREAMOutput(outfile)

                if deleteOutput:
                    os.remove(outfile)
        except TimeoutExpired:
            p.kill()
            errorOnExit = 3
    except KeyboardInterrupt:
        errorOnExit = 2
    finally:
        os.remove(infile)

    if errorOnExit == 1:
        print(stderr_data)
        raise DREAMException("DREAMi exited with a non-zero exit code: {}".format(p.returncode))
    elif errorOnExit == 2:
        raise DREAMException("DREAMi simulation was cancelled by the user.")
    elif errorOnExit == 3:
        raise DREAMException("DREAMi simulation was killed due to timeout.")
    else:
        return obj

locatedream()

