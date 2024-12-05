function h=precalc_intxn(h,p)

param = getparam(h);
intxn = calc_intxn(h,0:p);
param.intxn = intxn ; 
h = setparam(h,param);

