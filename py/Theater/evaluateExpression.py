# Helper module for evaluating user-provided Python code

def evaluate(code, output):
    """
    Evaluate the given Python code.

    :param str code: Python code to evaluate.
    :param DREAMOutput output: DREAMOutput object.
    """
    localvariables = {
        'do': output
    }
    exec(code, globals(), localvariables)


