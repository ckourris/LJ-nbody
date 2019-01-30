""" * Authors: C. Kourris and Ethan van Woerkom
    * This module implements the Main function and
    * program for an argon N-body simulation.
"""
import numpy as np
import sys
from Box import *
from Particle3D import Particle3D
from Utilities import *
from MDUtilities import *

def main(arg1, arg2):
    if len(sys.argv)!= 3:
        print("Wrong number of arguments.")
        print("Usage: " + sys.argv[0] + " <output file>")
        quit()
    else:
        param_name = sys.argv[1]
        out_name = sys.argv[2]

    parameters, outfile = get_arguments(param_name, out_name)

    Box = Box(parameters[0], parameters[2], parameters[1], parameters[3])

    Box.simulate()

    # Need the RSD etc methods here

main()
