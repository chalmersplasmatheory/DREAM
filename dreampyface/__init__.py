import sys
from pathlib import Path
p = str((Path(__file__).parent/'../build/dreampyface/cxx').resolve().absolute())
sys.path.append(p)
from .libdreampy import *
from .Simulation import Simulation


