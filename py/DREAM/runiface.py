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
    print("locatedream called")
    try:
        DREAMPATH = os.environ['DREAMPATH']

        if DREAMPATH[-1] == '/':
            DREAMPATH = DREAMPATH[:-1]

    except KeyError: pass

    if DREAMPATH is None:
        DREAMPATH = (pathlib.Path(__file__).parent / '..' / '..').resolve().absolute()
        print(DREAMPATH)
        if not os.path.isfile('{}/build/iface/dreami'.format(DREAMPATH)):
            #raise DREAMException("Unable to locate the DREAMi executable. Try to set the 'DREAMPATH' environment variable.")
            print("WARNING: Unable to locate the DREAMi executable. Try to set the 'DREAMPATH' environment variable.")

def runiface(settings, outfile=None, quiet=False, timeout=None, nthreads=None):
    """
    Run 'dreami' with the specified settings (which may be either
    a 'DREAMSettings' object or the name of a file containing the
    settings).

    :param settings: ``DREAMSettings`` object or name of file containing settings.
    :param outfile:  Name of file to write output to (default: 'output.h5')
    :param quiet:    If ``True``, suppresses output to stdout/stderr. (default: ``False``)
    :param timeout:  Time (in seconds) after which process should be killed. If ``None``, no time limit is imposed.
    :param nthreads: Number of OpenMP threads to allow at most (only relevant for the MKL linear solver).
    """
    global DREAMPATH
    task = DREAMTask(settings, outfile, quiet, timeout, DREAMPATH=DREAMPATH, nthreads=nthreads)
    task.run()
    while not task.hasFinished():
        pass
    return task.getResult()

def runiface_parallel(settings, outfiles,stdout_list=[],stderr_list=[], quiet=False, timeout=None, njobs=4, nthreads=None):
    global DREAMPATH
    queue = []
    active = []
    allTasks = []
    
    if len(settings) != len(outfiles):
        raise DREAMException("Lengths of settings and outfiles arrays are different!")
    
    #for _settings, outfile, stdout_name, stderr_name in zip(settings, outfiles,stdout_list,stderr_list):
    for i in range(len(settings)):
        so = '' if len(stdout_list)==0 else stdout_list[i]
        se = '' if len(stderr_list)==0 else stderr_list[i]
        task = DREAMTask(settings[i], outfiles[i], quiet, timeout, DREAMPATH, so, se, nthreads=nthreads)
           
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

    return [task.getResultObject() for task in allTasks]
    

locatedream()

