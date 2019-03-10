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
    parameters, outfile, cpp = get_arguments()

    # Create simulation Box. See design document for details.
    Simba = Box(parameters[0], parameters[2], parameters[1], parameters[3])
    position_list, timelist = Simba.simulate(outfile, parameters[5], parameters[4])

    # If you only want to load data to test the observable use the following
    # command and comment the above two. In that case specify the outfile:
    #outfile = ""
    #position_list = np.array(get_output(outfile, parameters[0]))

    # Define MSD and RDF parameters
    msd_start = 1 # Need something >= 1
    msd_end = 200 # Need something > msd_start and < n_steps
    rdf_bins = np.arange(0,5,0.1)
    rdf_start = 100 # Need > 0
    rdf_end = 200 # Need < n_steps

    print("Calculating the Mean Square Displacement function\n")
    MSD_arr = MSD(position_list, msd_start, msd_end, Simba.boxdim)
    #fig = plt.figure(figsize=(3, 6))
    write_output("MSD_output.txt", timelist[msd_start - 1:msd_end], MSD_arr)
    plt.plot(timelist[msd_start - 1:msd_end], MSD_arr)
    plt.title("Mean Square Displacement")
    plt.xlabel("Time")
    plt.ylabel("MSD(t)/sigma^2")
    #plt.show()
    plt.savefig('MSD.png')
    plt.close()

    print("Calculating the Radial Distribution function\n")
    rdf_arr, rdf_bins = RDF(position_list, rdf_start, rdf_end, rdf_bins, Simba.boxdim)
    write_output("RDF_output.txt", rdf_bins, rdf_arr)
    plt.plot(rdf_bins,rdf_arr)
    plt.title("Radial Distribution Function")
    plt.xlabel("Distance/sigma")
    plt.ylabel("g(r)")
    #plt.show()
    plt.savefig('RDF.png')

    print("Both plots saved in directory. Simulation has been successful")

if __name__ == '__main__':
    main()
