""" * Authors: C. Kourris and Ethan van Woerkom
    * This module implements various functions
    * for an argon N-body simulation.
"""
import numpy as np
import sys.argv

def get_arguments():
	""""The program is called with two filenames as input. The first one
	contains the simulation parameters. The second one is the name of
	the file to which the VMD output will be written.
	This function parses those inputs and returns a list with the
	contents of the first file
	and a string with the name of the second file."""
	if len(sys.argv) != 3:
		print("Wrong number of arguments, give two, e.g.:")
		print("Main.py parameters.txt vmdoutput.xyz")
		raise Exception('Wrong arguments')
	paramfilename = sys.argv[1]
	VMDfile = sys.argv[2]
	# Now get parameters individually:
	paramfile = open(paramfilename, 'r')
	lines = paramfile.readlines()
	N = int(lines[1])
	rho = float(lines[3])
	LJ_cutoff = float(lines[5])
	T = int(lines[7])
	dt = float(lines[]9)
	nsteps = int(lines[11])
	return param = [N, rho, LJ_cutoff, T, dt, nsteps], VMDfile


def LJ_Potential(vector, cutoff):
	""""This function calculates the Lennard-Jones potential for a
	particle separated by the vector in distance,
	using a given cutoff distance. It returns a scalar."""
    return scalarpotential


def LJ_Force(vector, cutoff):
	""" function calculates the Lennard-Jones force vector for a
	particle separated by the vector in distance,
	using a given cutoff distance. It returns an narray."""
    return force

def Total_PE(positions, cutoff):
	"""This function returns the total calculated potential energy
	for a set of system positions given by an [N, 3]-dimensional
	narray, using a given LJ cutoff."""
    return energy


def Total_KE(velocities):
	"""This function returns the total calculated kinetic energy
	for an [N,3] dimensional narray of system velocities."""
    squares = np.linalg.norm(velocities) # sum of squares of velocities
    mass = 1 # Setting argon mass to 1
    return 0.5*mass*squares


def RDF(pos, start, end, bins):
	"""Given a [T, N,3]-dimensional narray of system positions indexed
	by time, this will calculate the radial density function histogram
	averaged from time start to end exclusive using the given bins."""
    return radial_density_histogram


def MSD(pos, start, length):
	"""Given a [T, N,3]-dimensional narray of system positions indexed
	by time, this will calculate the mean square displacement for the
	system from time start to start+length exclusive
	relative to the given start time."""
    return mean_square_displacement
