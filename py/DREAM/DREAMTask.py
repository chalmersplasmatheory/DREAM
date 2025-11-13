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
    def __init__(self, settings, outfile=None, quiet=False, timeout=None, DREAMPATH=None, stdout_name='', stderr_name='', nthreads=None):
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
        #self.stderr_dat_byte=None
        self.DREAMPATH = DREAMPATH
        self.quiet = quiet
        self.timeout = timeout if timeout is not None else float('inf')
        self.nthreads = nthreads
        
        self.stdout_name =stdout_name   
        self.stderr_name =stderr_name  #new
        self.open_stdout = None
        self.open_stderr = None

    def run(self):
        self.startTime= time.time()
        if self.p != None: #if process is already created then we can safely od nothing
            return 
            
        if self.stdout_name != '':
                self.open_stdout=open(self.stdout_name,'w')
        else :  
                self.open_stdout=subprocess.PIPE

        if self.stderr_name != '':
                self.open_stderr=open(self.stderr_name,'w')
        else :  
                self.open_stderr=subprocess.PIPE

        env = None
        if self.nthreads is not None:
            env = os.environ.copy()
            env['OMP_NUM_THREADS'] = str(self.nthreads)
           
        if self.quiet: 
            self.p = subprocess.Popen([f'{self.DREAMPATH}/build/iface/dreami', self.infile], stderr=self.open_stderr, stdout=self.open_stdout, env=env)
        else:
            self.p = subprocess.Popen([f'{self.DREAMPATH}/build/iface/dreami', self.infile], stderr=self.open_stderr, env=env)

    def hasFinished(self, timeout=1):
        try:
            self.stderr_dat_byte=self.p.communicate(timeout=timeout)[1]
            if self.stderr_dat_byte != None:
                self.stderr_data = self.stderr_dat_byte.decode('utf-8')
           
            if self.stdout_name!='':
                self.open_stdout.close()
            if self.stderr_name!='':
                self.open_stderr.close()
                
            if self.p.returncode != 0:
                self.errorOnExit = 1
            else:  #gives error when closing file since outfiles do not correspond to the actual output directory files
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
