""" * Authors: C. Kourris and Ethan van Woerkom
    * This module implements the Box class and
    * methods for an argon N-body simulation.
"""
import Particle3D
import MDUtilities
from Utilities import LJ_Force
import numpy as np

class Box:
    """ CLASS VARIABLES:
    particles - List containing information about all particles.
    boxdim - Dimension of box.
    LJ_cutoff - Lennard-Jones cutoff distance.
    """
    def __init__(self, N, LJ_cutoff, rho, T):
        """
        Initialises simulation box with given parameters using
        function from MDUtilities.py to set particle positions and velocities.
        """
        # Initialise list of particles with zero position and velocity
        # and label equal to their number
        self.particles = [Particle3D.Particle3D(str(i)) for i in range(1,N+1)]
    
        self.LJ_cutoff = LJ_cutoff # Save LJ_cutoff distance.
    
        # Set particle positions, get box dimensions:
        self.boxdim = MDUtilities.set_initial_positions(rho, self.particles)[0]
        MDUtilities.set_initial_velocities(T, self.particles) # Set velocities
        return None

    
    def update_vel(self, forces, dt):
        """
        Conducts first order velocity update given
        an narray of forces on all particles, using v = v + F*d
        t"""
        for particle, force in zip(self.particles, forces):
            particle.leap_velocity(dt, force)
        return None
    
    def update_pos(self, forces, dt):
        """
        Conducts second order position update given an narray
        of forces on all particles, using x = x + v*dt + 0.5*F*dt^2.
        """
        for particle, force in zip(self.particles, forces):
            particle.leap_position(dt, force)
        return None


    def get_positions(self):
        """ Returns [N,3]-dim narray of positions of all particles. """
        return np.array([p.position for p in self.particles])


    def get_velocities(self):
        """ Returns [N,3]-dim narray of velocities of all particles."""
        return np.array([p.velocity for p in self.particles])


    def get_forces(self):
        """Returns [N,3]-dim narray of forces on all particles."""
        return particle_forces

    
    def VMD_string(self, time):
        """
        Produces a string in the VMD format giving the state of the
        current system, with Point = T. Uses the Particle3D.__str__ method.
        """
        return VMD_output_string

    def enforce_pbc(self):
        """
        Enforces period boundary conditions and moves any particle that
        has strayed outside the box back into the box according to pbc.
        """
        return None

    
    def simulate(self, outputfile, nsteps, dt):
        """
        Runs a Verlet n-body simulation on the initialised box for nsteps
        with timestep dt, and returns [nsteps,N,3]-dim position and velocity
        narrays and a nsteps-length time narray.
        """
        return positions, velocities, times


