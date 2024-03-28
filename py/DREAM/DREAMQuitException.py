# Implementation of a basic exception class
from .DREAMException import DREAMException

class DREAMQuitException(DREAMException):
    
    def __init__(self, message):
        super(Exception, self).__init__(message)


