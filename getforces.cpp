#include <math.h>
#include <iostream>

using namespace std;

double mod(double num, double div){
    double res = fmod(num, div);
    return (res >= 0) ? res : res+div;
}

void force(double* p1, double* p2, double* output, double boxdim, double cutoff){
    // First calculate MIC vector
//    cout << "p2[0] = " << p2[0] << endl;
    double sep[3];
    for(int i = 0; i < 3; i++) sep[i] = p2[i]-p1[i];
    for(int i = 0; i < 3; i++) sep[i] = mod(sep[i], boxdim);
    for(int i = 0; i < 3; i++) sep[i] = mod(sep[i]+boxdim/2, boxdim) - boxdim/2;
//    cout << "sep[0] = " << sep[0] << endl;
    // Calculate radius
    double r = sqrt(pow(sep[0], 2) + pow(sep[1], 2) + pow(sep[2], 2));
//    cout << "r = " << r << endl;
    // calculate force
    for(int i = 0; i < 3; i++){
        output[i] = (r<cutoff) ? 48*(pow(r,-14)-0.5*pow(r,-8))*sep[i] : 0;
    }

    return;
}

void addvector(double* invec, double* outvec, bool negative){
    for(int i = 0; i < 3; i++) outvec[i] += invec[i]* (negative ? -1 : 1);
    return;
}

void getforces(double* in_array, double* out_array, int N, double boxdim, double cutoff){
    //Set out ndarray to zero
    for(int i = 0; i<3*N; i++) out_array[i]=0;

    // Setup temporary force vector;
    double forcevar[3];

    // Calculate all forces by iterating over unique i, j pairs
    for(int i = 0; i<N; i++){
        for(int j = 0; j < i; j++){
            force(in_array+i*3, in_array+j*3, forcevar, boxdim, cutoff);
//          cout << forcevar[0] << ' ' << forcevar[1] << endl;
            addvector(forcevar, out_array+j*3, 0);
            addvector(forcevar, out_array+i*3, 1);
        }
    }

    return;
}
