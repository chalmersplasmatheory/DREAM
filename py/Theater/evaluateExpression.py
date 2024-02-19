# Helper module for evaluating user-provided Python code

import matplotlib.pyplot as plt
import numpy as np


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


