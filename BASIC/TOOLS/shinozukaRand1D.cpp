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

    // V = shinozukaRand1D(order,nx,c,z,phi,k,x);

    mwSize order, nx;
    order = (mwSize) *mxGetPr( prhs[ 0 ] );
    nx = (mwSize) *mxGetPr( prhs[ 1 ] );

    double *c;
    c = (double*) mxGetPr( prhs[ 2 ] );
    double *z;
    z = (double*) mxGetPr( prhs[ 3 ] );
    double *phi;
    phi = (double*) mxGetPr( prhs[ 4 ] );
    double *k;
    k = (double*) mxGetPr( prhs[ 5 ] );

    void *p;
    p = (void *) mxGetPr( prhs[ 6 ] );
    Map< Array < double, Dynamic, 1 > > x( (double *) p, nx, 1 );

    const mwSize dims[]={nx};
    plhs[0] = mxCreateNumericArray( 1, dims, mxDOUBLE_CLASS, mxREAL );

    p = (void *) mxGetData( plhs[0] );
    Map< Array < double, Dynamic, 1 > > V( (double *) p, nx, 1 );

    V.setZero();

    double PI = 3.141592653589793;

    ArrayXXd Vi(nx,1);

    for(int i=0; i<order; i++)
    {
//         cout << k1[i] << k2[i] << endl;
        Vi = sqrt(2) * c[i] * z[i] * cos(phi[i] + k[i]*x);
        V += Vi;
    }

}
