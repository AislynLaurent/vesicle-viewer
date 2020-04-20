# Bring in functions
from math import *
import numpy as np

def safe_function_dict():
    # Declasafe functions for eval
    safe_function_dict = {
        'np.piecewise':np.piecewise,
        'np.e':np.e,
        'np.floor':np.floor,
        'np.ceil':np.ceil,
        'np.log':np.log,
        'exp':exp,
        'np.log10':np.log10,
        'np.pi':np.pi,
        'pow':pow,
        'np.sqrt':np.sqrt,
        'abs':abs,
    }

    return safe_function_dict