/*
 * fastprodscal.c
 * multiply the vectors of a SEPMATRIX
 * To compile the function use :
 *
 * For linux
 * mex -largeArrayDims fastprodscal.c -lmwblas
 *
 * For Windows
 * mex('fastprodscal.c',[matlabroot '\extern\lib\win32\lcc\libmwblas.lib'],'-largeArrayDims')
 *
 */

#include "mex.h"
#include "blas.h"
#include "matrix.h"

/* The C MEX-file gateway routine: */
void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]){
    int IndexOne=0;
    const mxArray *pu = prhs[0]; 
    const mxArray *puf = mxGetField(pu,IndexOne,"F");
    const mxArray *pua = mxGetField(pu,IndexOne,"alpha");
    
    size_t um = mxGetM(puf);
    size_t udim = mxGetN(puf);
    bool uiscell = mxIsCell(puf);
    
    const mxArray *pv = prhs[1];
    const mxArray *pvf = mxGetField(pv,IndexOne,"F");
    const mxArray *pva = mxGetField(pv,IndexOne,"alpha");
    
    size_t vm = mxGetM(pvf);
    size_t vdim = mxGetN(pvf);
    bool viscell = mxIsCell(pvf);
    
    const mxArray *ufid;
    const mxArray *vfjd;
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
    
    
    /* Check the number of dimensions */
    if (udim != vdim){
        mexErrMsgTxt("MULTIPLYF: the two inputs must have the same number of dimensions.");
    }
    
    
    double pij,PROD;
    double *uaid,*vajd,*UiVj;
    mwSignedIndex p,one;
    one=1;
    uaid = mxGetPr(pua);
    vajd = mxGetPr(pva);
    mxArray *rhs[2];
    PROD=0;
    for(i=0; i<um; ++i){
        for(j=0; j<vm; ++j){
            pij=uaid[i]*vajd[j];
            for(d=0; d<udim; ++d){
                ufid = mxGetCell(puf, d*um+i);
                vfjd = mxGetCell(pvf, d*vm+j);
                
                rhs[0] = ufid;
                rhs[1] = vfjd;
                mexCallMATLAB(1, &wfc, 2, rhs, "fastprodscal");
                UiVj = mxGetPr(wfc);
                pij  = pij * UiVj[0];
            }
            PROD = PROD + pij;
        }
    }
    /* OUTPUT */
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    *mxGetPr(plhs[0])=PROD;
}






