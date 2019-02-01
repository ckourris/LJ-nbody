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
    # Read parameter and output file names from sys.argv
    parameters, outfile = get_arguments()

    Box = Box(parameters[0], parameters[2], parameters[1], parameters[3])

    Box.simulate()

    # Need the RSD etc methods here

if __name__ == '__main__':
    main()
