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

def main():
    # Read parameter and output file names from sys.argv
    parameters, outfile = get_arguments()

    Simba = Box(parameters[0], parameters[2], parameters[1], parameters[3])

    Simba.simulate(outfile, parameters[5], parameters[4])

    with open(outfile, 'r') as f:
         position_list = np.array([line.strip().split() for line in f if not \
         (line.startswith(('P')) or line.startswith(('108')))])

    MSD_arr = MSD(position_list, 1,100, Simba.boxdim)

if __name__ == '__main__':
    main()
