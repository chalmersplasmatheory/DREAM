

class DREAMErrorOutput:

    def __init__(self, exception, output = None, errorCode = -1):
        """
        Object used to hold a result of invalid DREAMTask execution.
        """
        self.exception = exception
        self.output = output
        self.errorCode = errorCode