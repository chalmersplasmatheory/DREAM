# Implementation of a basic exception class

class EquationException(Exception):
    
    def __init__(self, message):
        super(Exception, self).__init__(message)


