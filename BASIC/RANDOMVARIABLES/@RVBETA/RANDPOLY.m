function p = RANDPOLY(r)

param=getparam(r);
p=POLYJACOBI(param.a,param.b);
p=setnumber(p,getnumber(r));
