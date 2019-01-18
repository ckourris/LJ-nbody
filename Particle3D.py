""" * Authors: C. Kourris and Ethan van Woerkom
    * This module implements the Particle3D class and
    * methods for an argon N-body simulation.
"""
import numpy as np

class Particle3D(object):
    """
    Class to describe 3D particles.

    Properties:
    label - a string describing the type of the particle
    position - float valued (numpy) array of position eg: [x,y,z]
    velocity - float valued (numpy) array of velocity eg: [vx,vy,vz]
    mass(float) - particle's mass

    Methods:
    * formatted output
    * kinetic energy
    * first-order velocity update
    * first- and second order position updates
    * A static method to create a particle from a given file entry
    * A static method that returns the relative separation vector
      between two particle objects.
    """

    def __init__(self, label, pos, vel, mass):
        """
        Initialise a Particle3D instance

        :param label: label as string
        :param pos: position as float valued array
        :param vel: velocity as float valued array
        :param mass: mass as float
        """

        self.position = pos
        self.velocity = vel
        self.mass = mass
        self.label = label

    def __str__(self):
        """
        Define output format.
        For particle p=(p_label 0.0 1.0 2.0 0.5 0.5 0.5 1) this will print as
        "p_label x = 0.0, y = 1.0, z = 2.0"
        """

        return self.label + " x = " + str(self.position[0]) + ", y = " + str(self.position[1]) + ", z = " + str(self.position[2])

    def kinetic_energy(self):
        """
        Return kinetic energy as
        1/2*mass*vel^2
        """

        speed = np.linalg.norm(self.velocity)
        return 0.5*self.mass*speed**2

    def leap_velocity(self, dt, force):
        """
        First-order velocity update,
        v(t+dt) = v(t) + dt*F(t)/m

        :param dt: timestep as float
        :param force: force on particle as float
        """

        self.velocity = self.velocity + dt*force/self.mass

    def leap_pos1st(self, dt):
        """
        First-order position update,
        x(t+dt) = x(t) + dt*v(t)

        :param dt: timestep as float
        """

        self.position = self.position + dt*self.velocity

    def leap_pos2nd(self, dt, force):
        """
        Second-order position update,
        x(t+dt) = x(t) + dt*v(t) + 1/2*dt^2*F(t)

        :param dt: timestep as float
        :param force: current force as float
        """

        self.position = self.position + dt*self.velocity + 0.5*dt**2*force/self.mass

    @staticmethod
    def from_file(file_handle):
        """
        A static method to create a particle from a file entry. The method takes a line
        and splits it into the arguments of an object. The expected form of the line is:
        "label pos_x pos_y pos_z vel_x vel_y vel_z mass"

        :param: An open file that can be accessed and read
        :return: An instance of an object (a particle object)
        """

        line = file_handle.readline()
        args = line.split(" ")
        label = args[0]
        position = np.array(args[1:4],float)
        velocity = np.array(args[4:7],float)
        mass = float(args[7])
        return Particle3D(label,position,velocity,mass)

    @staticmethod
    def relative_sep(p1, p2):
        """
        Return the relative separatin position vector directed from p1 to p2

        :param p1: A particle object
        :param p2: A particle object
        :return: Separation vector as numpy array
        """

        return p1.position - p2.position
