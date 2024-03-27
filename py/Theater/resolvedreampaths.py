# Resolve, and if necessary, add paths needed for the
# various required DREAM Python libraries.

import pathlib
import sys


cwd = pathlib.Path(__file__).resolve().parent.resolve()


# Add missing paths if necessary
try: import DREAM
except ModuleNotFoundError:
    sys.path.append(str(cwd.parent.resolve()))

try: import dreampyface
except ModuleNotFoundError:
    sys.path.append(str(cwd.parent.parent))

try: import libdreampy
except ModuleNotFoundError:
    sys.path.append(str(cwd.parent.parent / "build/dreampyface/cxx"))

