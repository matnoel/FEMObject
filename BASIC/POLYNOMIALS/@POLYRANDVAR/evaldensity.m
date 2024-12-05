function px=evaldensity(h,x)
param=getparam(h);
px=pdf(param.r,x);
