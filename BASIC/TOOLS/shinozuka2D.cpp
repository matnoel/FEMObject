#include "mex.h"
#include "matrix.h"
#include <math.h>
#include <iostream>
#include <unistd.h>

#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    //unsigned nb_core = (unsigned) *mxGetPr( prhs[ 10 ] );
    //setNbThreads(nb_core);
    //int nb_threads = Eigen::nbThreads( );
    //cout << "Number of threads (sub_Lrond_mex) = " << nb_threads << endl ;

//     int numCPU = sysconf(_SC_NPROCESSORS_ONLN);
//     setNbThreads(numCPU);

    //cout << "Number of cores (sub_Lrond_post_in_mex) = " << numCPU << endl ;

    // Transfert data between matlab and C++

    // V = shinozuka2D(nu,nx,c1,c2,z,phi,k1,k2,x1,x2);

    mwSize nu, nx;
    nu = (mwSize) *mxGetPr( prhs[ 0 ] );
    nx = (mwSize) *mxGetPr( prhs[ 1 ] );

    double *c1;
    c1 = (double*) mxGetPr( prhs[ 2 ] );
    double *c2;
    c2 = (double*) mxGetPr( prhs[ 3 ] );
    double *z;
    z = (double*) mxGetPr( prhs[ 4 ] );
    double *phi;
    phi = (double*) mxGetPr( prhs[ 5 ] );
    double *k1;
    k1 = (double*) mxGetPr( prhs[ 6 ] );
    double *k2;
    k2 = (double*) mxGetPr( prhs[ 7 ] );

    void *p;
    p = (void *) mxGetPr( prhs[ 8 ] );
    Map< Array < double, Dynamic, 1 > > x1( (double *) p, nx, 1 );
    p = (void *) mxGetPr( prhs[ 9 ] );
    Map< Array < double, Dynamic, 1 > > x2( (double *) p, nx, 1 );

    const mwSize dims[]={nx};
    plhs[0] = mxCreateNumericArray( 1, dims, mxDOUBLE_CLASS, mxREAL );

    p = (void *) mxGetData( plhs[0] );
    Map< Array < double, Dynamic, 1 > > V( (double *) p, nx, 1 );

    V.setZero();

    double PI = 3.141592653589793;

    ArrayXXd Vij(nx,1);

    for(int i=0; i<nu; i++)
    {
        for(int j=0; j<nu; j++)
        {
//             cout << k1[i] << k2[j] << endl;
            Vij = sqrt(2) * c1[i]*c2[j] * z[i + j*nu] * cos(phi[i + j*nu] + k1[i]*x1 + k2[j]*x2);
            V += Vij;
        }
    }

}
