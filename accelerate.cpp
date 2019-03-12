/*  * Authors: C. Kourris and Ethan van Woerkom
    * This module implements the forces and energies c++
    * accelerated functions and is imported as a library
    * into the Python code.
    * Note: Compilation might be necessary. See README file
*/

#include <math.h>

using namespace std;

double mod(double num, double div){
    // Implements a mod function for floating point numbers that always maps
    // to the range [0,div).

    double res = fmod(num, div);
    return (res >= 0) ? res : res+div;
}


void force(double* p1, double* p2, double* output, double boxdim, double cutoff){
    /* Calculates the force between two particles and writes result to output
     * Param:
     * p1, p2: Arrays of length 3 representing the position vectors
     * output: Array of length 3 where output vector is to be written.
     * boxdim: double representing dimensions of PBC box.
     * cutoff: double representing the LJ cutoff distance
     */

    // Calculates the MIC separation vector and stores it in sep.
    double sep[3];
    for(int i = 0; i < 3; i++) sep[i] = p2[i]-p1[i];
    for(int i = 0; i < 3; i++) sep[i] = mod(sep[i], boxdim);
    for(int i = 0; i < 3; i++) sep[i] = mod(sep[i]+boxdim/2, boxdim) - boxdim/2;

    // Calculate distance between particles.
    double r = sqrt(pow(sep[0], 2) + pow(sep[1], 2) + pow(sep[2], 2));

    // calculate force
    for(int i = 0; i < 3; i++){
        output[i] = (r<cutoff) ? 48*(pow(r,-14)-0.5*pow(r,-8))*sep[i] : 0;
    }

    return;
}


double potential(double* p1, double* p2, double boxdim, double cutoff){
  /* Calculates the potential energy between two particles and returns the value.
   * Param:
   * p1, p2: Arrays of length 3 representing the position vectors
   * boxdim: double representing dimensions of PBC box.
   * cutoff: double representing the LJ cutoff distance
   */

    // Calculates the MIC separation vector and stores it in sep.
    double sep[3], output = 0;
    for(int i = 0; i < 3; i++) sep[i] = p2[i]-p1[i];
    for(int i = 0; i < 3; i++) sep[i] = mod(sep[i], boxdim);
    for(int i = 0; i < 3; i++) sep[i] = mod(sep[i]+boxdim/2, boxdim) - boxdim/2;

    // Calculate distance between particles.
    double r = sqrt(pow(sep[0], 2) + pow(sep[1], 2) + pow(sep[2], 2));

    // Calculate potential energy between particles.
    output = (r<cutoff) ? 4*(pow(r,-12)-pow(r,-6))-4*(pow(cutoff,-12)-pow(cutoff,-6)) : 0;

    return output;
}


double KE(double* v){ //, double* out){
  /* Calculates the Kinetic energy corresponding to velocity vector stored
   * in the array v of length 3, K = 1/2*m*v^2, and return the value.
   */

  double temp_k = 0;

  for (int i = 0; i < 3; i++){
    temp_k += 0.5*(double)pow(v[i], 2);
  }
  return temp_k;
}


void addvector(double* invec, double* outvec, bool negative){
  /* Adds the vector invec of length 3 to outvec of length 3. Both are arrays.
   * The parameter negative determines the sign of the addition. If True invec
   * is subtracted from outvec.
   */

    for(int i = 0; i < 3; i++) outvec[i] += invec[i]* (negative ? -1 : 1);

    return;
}


void getforces(double* in_array, double* out_array, int N, double boxdim, double cutoff){
  /* Calculates the forces between N particles and writes the results to out array.
   * Param:
   * in_array: Array of dimensions (N,3) that holds the positions of N particles.
   * out_array: Array of dimensions (N,3) in which the forces of N particles are
   *            are written.
   * N: Number of particles.
   * boxdim: double representing dimensions of PBC box.
   * cutoff: double representing the LJ cutoff distance
   */

    //Set output numpy array to zero
    for(int i = 0; i<3*N; i++) out_array[i]=0;
    // Setup temporary force vector;
    double forcevar[3];

    // Calculate all forces by iterating over unique i, j pairs
    for(int i = 0; i<N; i++){
        for(int j = 0; j < i; j++){
            force(in_array+i*3, in_array+j*3, forcevar, boxdim, cutoff);
            addvector(forcevar, out_array+j*3, 0);
            addvector(forcevar, out_array+i*3, 1); // F_reaction = -F_action
        }
    }

    return;
}


void getenergies(double* pos_array, double* v_array, double* out_array,
                  int N, double boxdim, double cutoff){
    /* Calculates the potential, kinetic and total energy between N particles
     * and writes the results to out array.
     * Param:
     * pos_array: Array of dimensions (N,3) that holds the positions of N particles.
     * v_array: Array of dimensions (N,3) that holds the velocities of N particles.
     * out_array: Array of length 3 to which Potential, Kinetic and Total energy are
     *            written, in that order.
     * N: Number of particles.
     * boxdim: double representing dimensions of PBC box.
     * cutoff: double representing the LJ cutoff distance
     */

    double kinetic = 0;
    double poten = 0;

    // Calculate Potential energy for all pairs of particles
    for (int i = 0; i<N; i++){
      for (int j = 0; j < i; j++){
        poten += potential(pos_array+i*3, pos_array+j*3, boxdim, cutoff);
      }
    }

    // Calculate Kinetic energy for all particles.
    for (int i = 0; i<N; i++){
        kinetic += KE(v_array+i*3);
    }

    // Store results in return array.
    out_array[0] = poten;
    out_array[1] = kinetic;
    out_array[2] = poten+kinetic;

    return;
}
