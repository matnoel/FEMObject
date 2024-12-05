function rv=RANDVAR(h)
param=getparam(h);
rv=RVGAMMA(param.a,1,0);
rv = setnumber(rv,getnumber(h));