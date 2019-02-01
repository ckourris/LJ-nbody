""" * Authors: C. Kourris and Ethan van Woerkom
    * This module implements the Box class and
    * methods for an argon N-body simulation.
"""
from Particle3D import Particle3D
import MDUtilities
from Utilities import LJ_Force
import numpy as np
import time
import sys

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
        self.particles = [Particle3D(str(i)) for i in range(1,N+1)]

        self.LJ_cutoff = LJ_cutoff # Save LJ_cutoff distance.

        # Set particle positions, get box dimensions:
        self.boxdim = MDUtilities.set_initial_positions(rho, self.particles)[0]
        MDUtilities.set_initial_velocities(T, self.particles) # Set velocities
        return None


    def update_vel(self, forces, dt):
        """
        Conducts first order velocity update given
        an narray of forces on all particles, using v = v + F*dt
        """
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
        N = len(self.particles)
        particle_forces = np.zeros( (N,3) )

        # Iterate over all i<j, then calculate
        # force for each i, j combination
        for i in range(N):
            for j in range(i):
                # Get force of particle i on j, respecting pbc and mic.
                sep = Particle3D.pbc_sep(self.particles[i], self.particles[j], self.boxdim)
                force = LJ_Force(sep, self.LJ_cutoff)
                particle_forces[j] += force
                particle_forces[i] += -force # Using Newtons 3rd law

        return particle_forces

    def VMD_string(self, time):
        """
        Produces a string in the VMD format giving the state of the
        current system, with Point = T. Uses the Particle3D.__str__ method.
        """
        # First add preamble:
        # N_data
        # Point = time
        # ...
        VMD_output_string = '%i\nPoint = %i\n' % (len(self.particles), time)
        # Concatenate labels and positions for each particle
        for p in self.particles:
            VMD_output_string += str(p) + '\n'
        return VMD_output_string

    def enforce_pbc(self):
        """
        Enforces period boundary conditions and moves any particle that
        has strayed outside the box back into the box according to pbc.
        """
        for particle in self.particles:
            particle.position = np.mod(particle.position,self.boxdim)
        return None


    def simulate(self, outputfile, nsteps, dt):
        """
        Runs a Verlet n-body simulation on the initialised box for nsteps
        with timestep dt, and returns [nsteps,N,3]-dim position and velocity
        narrays and a nsteps-length time narray.
        """
        starttime = time.process_time()
        timelist = []; VMD_list = []; positions = []; velocities = [];
        forces = self.get_forces()
        for t in range(nsteps):
            positions.append(self.get_positions())
            self.enforce_pbc()
            velocities.append(self.get_velocities())
            timelist.append(t*dt)
            VMD_list.append(self.VMD_string(t))

            # Updates positions, velocities etc
            self.update_pos(forces, dt)
            temp_forces = forces
            forces = self.get_forces()
            self.update_vel(0.5*(temp_forces + forces), dt)

        # Output VMD data to file
        vmdstring = ''.join(VMD_list)
        with open(outputfile, 'w') as out:
            out.write(vmdstring)
            print('Succesful VMD Data write to '+outputfile)

        runtime = time.process_time() - starttime
        print('Simulate method ran for %f seconds'%runtime)
        return np.array(positions), np.array(velocities), np.array(timelist)
