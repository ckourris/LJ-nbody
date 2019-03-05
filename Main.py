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
import matplotlib.pyplot as plt
import time

def main():
    # Read parameter and output file names from sys.argv
    parameters, outfile = get_arguments()

    Simba = Box(parameters[0], parameters[2], parameters[1], parameters[3])
    Simba.simulate(outfile, parameters[5], parameters[4])

    position_list = np.array(get_output(outfile, parameters[0]))

    #MSD_arr,times = MSD(position_list, 1,500, Simba.boxdim)
    #plt.plot(times, MSD_arr)
    #plt.show()

    rdf_arr = RDF(position_list, 200, 300, np.arange(0,5,0.1), Simba.boxdim)
    rdf = Bin_particles(position_list[10], np.arange(0,5,0.1), Simba.boxdim)
    plt.plot(np.arange(0.1,5.0,0.1),rdf_arr)
    plt.show()

if __name__ == '__main__':
    main()
