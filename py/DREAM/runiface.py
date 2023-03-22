# Simple wrapper for running 'dreami'

import os
import pathlib
import subprocess
import tempfile

from . DREAMException import DREAMException
from . DREAMOutput import DREAMOutput
from . DREAMSettings import DREAMSettings
from . DREAMTask import DREAMTask
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
    task = DREAMTask(settings, outfile, quiet, timeout, DREAMPATH )
    task.run()
    while not task.hasFinished():
        pass
    return task.getResult()

def runiface_parallel(settings, outfiles, quiet=False, timeout=None, njobs=4):
    global DREAMPATH
    queue = []
    active = []
    allTasks = []
    if len(settings) != len(outfiles):
        raise DREAMException("Lengths of settings and outfiles arrays are different!")
    
    for _settings, outfile in zip(settings, outfiles):
        task = DREAMTask(_settings, outfile, quiet, timeout, DREAMPATH )
        queue.append(task)
        allTasks.append(task)
    
    while len(queue) > 0 or len(active) > 0:
        for task in active:
            if task.hasFinished(0.1):
                active.remove(task)
        
        while len(queue) > 0 and len(active) < njobs:
            task = queue.pop(0)
            task.run()
            active.append(task)

    return [task.getResult() for task in allTasks]
    

locatedream()

