function rv=RANDVAR(h)
param=getparam(h);
rv=RVBETA(param.a,param.b,-1,1);
rv = setnumber(rv,getnumber(h));