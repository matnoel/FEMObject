/*
 * redux.c
 *
 *
 * function M = redux(a,Z)
 * M = a(1)*Z{1};
 * for i=2:A.m
 *     M = M + a(i)*Z{i};
 * end
 * return 
 *
 *******************************************************
 * WARNING: IF Z contains sparse matrices,
 * redux supposes that they all have the same structure!
 * It is not a bug, it is a feature.
 *******************************************************
 *
 * To compile the function use :
 *
 * For linux
 * mex -largeArrayDims redux.c -lmwblas
 *
 * For Windows
 * mex('redux.c',[matlabroot '\extern\lib\win32\lcc\libmwblas.lib'],'-largeArrayDims')
 *
 *
 */

#include "mex.h"
#include "blas.h"
#include "matrix.h"

/* The C MEX-file gateway routine: */
void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]){
    
    const mxArray *pa, *pz, *zi;
    mxArray *pm, *pnelem;
    double *a,*z,*m;
    bool aIsDouble, zIsCell, zIsDouble, zIsSparse;
    int sza,szz;
    mwSize nelem;
    ptrdiff_t incx, incy; /*using ptrdiff_t instead of int avoid many many problems*/
    int i;
    
    if(nlhs != 1){
        mexErrMsgTxt("One output argument required.");
    }
    if(nrhs !=2){
        mexErrMsgTxt("Two input arguments required.");
    }
    
    pa=prhs[0];
    aIsDouble=mxIsDouble(pa);
    
    pz=prhs[1];
    zIsCell=mxIsCell(pz);
    
    if(!aIsDouble || !zIsCell){
        mexErrMsgTxt("Usage: redux(a,Z), a is double, Z is cell");
    }
    a=mxGetPr(pa);
    sza=mxGetNumberOfElements(pa);
    
    zi=mxGetCell(pz,0);
    zIsDouble=mxIsDouble(zi);
    szz=mxGetNumberOfElements(pz);
    if(szz!=sza){
        mexErrMsgTxt("Usage: redux(a,Z), a and Z must have the same number of elements");
    }
    if(!zIsDouble){
        mexErrMsgTxt("Usage: redux(a,Z), Z is cell of full or sparse matrices");
    }
    
    zIsSparse=mxIsSparse(zi);
    if (zIsSparse){
        nelem=*(mxGetJc(zi) + mxGetN(zi)); /*number of non zero elements*/
    }else{
        nelem=mxGetNumberOfElements(zi);
    }

    pm=mxDuplicateArray(zi); /* M = Z{1} */
    m=mxGetPr(pm);

    incx=1;
    incy=1;
    dscal(&nelem, a, m, &incx); /* M = a(1)*M */

    for (i=1;i<sza;++i){
        zi=mxGetCell(pz,i);
        z=mxGetPr(zi);
        daxpy(&nelem,a+i,z,&incx,m,&incy); /* M = a(i)*Z{i} + M */
    }
    plhs[0]=pm;
}
