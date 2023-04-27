import os
import pathlib
import subprocess
import tempfile
import time

from . DREAMException import DREAMException
from . DREAMOutput import DREAMOutput
from . DREAMSettings import DREAMSettings
from . DREAMErrorOutput import DREAMErrorOutput
from subprocess import TimeoutExpired

class DREAMTask:
    def __init__(self, settings, outfile=None, quiet=False, timeout=None, DREAMPATH=None):
        self.deleteOutput = False
        if outfile is None:
            self.deleteOutput = True
            self.outfile = next(tempfile._get_candidate_names())+'.h5'
        else:
            self.outfile = outfile

        self.infile = None
        if isinstance(settings, DREAMSettings):
            self.infile = next(tempfile._get_candidate_names())+'.h5'
            settings.output.setFilename(self.outfile)
            settings.save(self.infile)
        else:
            self.infile = settings

        self.settings = settings
        self.errorOnExit = 0
        self.p = None
        self.obj = None
        self.stderr_data = None
        self.DREAMPATH = DREAMPATH
        self.quiet = quiet
        self.timeout = timeout if timeout is not None else float('inf')

    def run(self):
        self.startTime= time.time()
        if self.p != None: #if process is already created then we can safely od nothing
            return 
        if self.quiet:
            self.p = subprocess.Popen(['{}/build/iface/dreami'.format(self.DREAMPATH), self.infile], stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        else:
            self.p = subprocess.Popen(['{}/build/iface/dreami'.format(self.DREAMPATH), self.infile], stderr=subprocess.PIPE)

    def hasFinished(self, timeout=1):
        try:
            self.stderr_data = self.p.communicate(timeout=timeout)[1].decode('utf-8')

            if self.p.returncode != 0:
                self.errorOnExit = 1
            else:
                self.obj = DREAMOutput(self.outfile)

                if self.deleteOutput:
                    os.remove(self.outfile)
        except TimeoutExpired: 
            # In this case it is expected situation. Process has not completed during specified timeout.
            # We don't need to preserve any state, because documentation says that no stream data will be lost.

            # Now we implement timeout manually, by checking start time with current time
            if time.time() - self.startTime > self.timeout:
                self.p.kill()
                self.errorOnExit = 3
                return True
            return False
        except KeyboardInterrupt:
            self.errorOnExit = 2

        return True

    def getResult(self):
        os.remove(self.infile)

        if self.errorOnExit == 1:
            print(self.stderr_data)
            raise DREAMException("DREAMi exited with a non-zero exit code: {}".format(self.p.returncode))
        elif self.errorOnExit == 2:
            raise DREAMException("DREAMi simulation was cancelled by the user.")
        elif self.errorOnExit == 3:
            raise DREAMException("DREAMi simulation was killed due to timeout.")
        else:
            return self.obj

    def getResultObject(self):
        # Handling errors in this way is easiest in terms of backwards compatibility, but
        # quite inefficient. We can change it later, if we face performance issues.
        try:
            return self.getResult()
        except DREAMException as error:
            return DREAMErrorOutput(error, self.stderr_data, self.errorOnExit)