# Implementation of a basic exception class

class OutputException(Exception):
    
    def __init__(self, message):
        super(Exception, self).__init__(message)


