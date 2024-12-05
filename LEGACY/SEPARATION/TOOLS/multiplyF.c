/*
 * multiplyF.c
 * multiply the vectors of a SEPMATRIX
 * To compile the function use :
 *
 * For linux
 * mex -largeArrayDims multiplyF.c -lmwblas
 *
 * For Windows
 * mex('multiplyF.c',[matlabroot '\extern\lib\win32\lcc\libmwblas.lib'],'-largeArrayDims')
 *
 */

#include "mex.h"
#include "blas.h"
#include "matrix.h"

/* The C MEX-file gateway routine: */
void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]){
    const mxArray *puf = prhs[0]; 
    size_t um = mxGetM(puf);
    size_t udim = mxGetN(puf);
    bool uiscell = mxIsCell(puf);
    const mxArray *pvf = prhs[1];
    size_t vm = mxGetM(pvf);
    size_t vdim = mxGetN(pvf);
    bool viscell = mxIsCell(pvf);
    const mxArray *ufi;
    const mxArray *vfj;
    mxArray *wfc;

      int c=0;
    
    int d,i,j;
    
  
    
    /* Check for proper number of arguments */
    if (nrhs != 2) {
        mexErrMsgTxt("Two input arguments required.");
    } else if (nlhs != 1) {
        mexErrMsgTxt("One output argument required.");
    }
    
    /* INPUTS */
    
    
    /* Check the nature of the input */
    if (!uiscell || !viscell ){
        mexErrMsgTxt("MULTIPLYF requires cells as inputs.");
    }
    
    /* Check the number of dimensions */
    if (udim != vdim){
        mexErrMsgTxt("MULTIPLYF: the two inputs must have the same number of dimensions.");
    }
    
    /* OUTPUT */
    plhs[0] = mxCreateCellMatrix(um*vm, udim);
    
    /* PROCESSING */
    for(d=0; d<udim; ++d){
        if( mxIsSparse(mxGetCell(puf,d*um)) || mxIsSparse(mxGetCell(pvf,d*vm)) ){
            mxArray *rhs[2];
            for(i=0; i<um; ++i){
                ufi = mxGetCell(puf, d*um+i);
                for(j=0; j<vm;++j){
                    vfj = mxGetCell(pvf, d*vm+j);

                    rhs[0] = ufi;
                    rhs[1] = vfj;
                    mexCallMATLAB(1, &wfc, 2, rhs, "mtimes");
                    mxSetCell(plhs[0], c, wfc);
                    c+=1; 
                } 
            }
        }else { /* u & v is full */
            double *ui, *vj, *wc;
            mwSignedIndex m, n, p;
            char *chn = "N";
            double one = 1.0, zero = 0.0;
            for(i=0; i<um; ++i){
                ufi = mxGetCell(puf, d*um+i);
                ui = mxGetPr(ufi);

                m = (mwSignedIndex)mxGetM(ufi);
                p = (mwSignedIndex)mxGetN(ufi);

                for(j=0; j<vm;++j){
                    vfj = mxGetCell(pvf, d*vm+j);
                    vj = mxGetPr(vfj);
                    n = (mwSignedIndex)mxGetN(vfj);

                    if (p != mxGetM(vfj)) {
                        mexErrMsgIdAndTxt("MATLAB::multiplyF:matchdims",
                                "Inner dimensions of matrix multiply do not match.");
                    }

                    /* create output matrix */
                    wfc = mxCreateDoubleMatrix(m, n, mxREAL);
                    wc = mxGetPr(wfc);

                    /* Pass arguments to Fortran by reference */
                    dgemm(chn, chn, &m, &n, &p, &one, ui, &m, vj, &p, &zero, wc, &m);

                    mxSetCell(plhs[0], c, wfc);
                    c+=1;
                }
            }
        }
    }
}


