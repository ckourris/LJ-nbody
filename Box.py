""" * Authors: C. Kourris and Ethan van Woerkom
    * This module implements the Box class and
    * methods for an argon N-body simulation.
"""

class Box:
    particles # List containing information about all particles.
    boxdim #Dimension of box.
    LJ_cutoff #Lennard-Jones cutoff distance.

    def __init__(self, N, LJ_cutoff, rho, T):
""" Initialises simulation box with given parameters using
    function from MDUtilities.py to set particle positions and velocities. """
    return None

    
    def update_vel(self, forces, dt):
""" Conducts first order velocity update given
    an narray of forces on all particles, using v = v + F*dt"""
        return None

    
    def update_pos(self, forces, dt):
""" Conducts second order position update given an narray
    of forces on all particles, using x = x + v*dt + 0.5*F*dt^2."""
        return None


    def get_positions(self):
""" Returns [N,3]-dim narray of positions of all particles. """
        return particles_positions


    def get_velocities(self):
""" Returns [N,3]-dim narray of velocities of all particles."""
        return particles_velocities


    def get_forces(self):
    """Returns [N,3]-dim narray of forces on all particles."""
        return particle_forces

    
    def VMD_string(self, time):
""" Produces a string in the VMD format giving the state of the
    current system, with Point = T. Uses the Particle3D.__str__ method."""
        return VMD_output_string

    def enforce_pbc(self):
""" Enforces period boundary conditions and moves any particle that
    has strayed outside the box back into the box according to pbc."""
        return None

    
    def simulate(self, outputfile, nsteps, dt):
""" Runs a Verlet n-body simulation on the initialised box for nsteps
    with timestep dt, and returns [nsteps,N,3]-dim position and velocity
    narrays and a nsteps-length time narray."""
        return positions, velocities, times


