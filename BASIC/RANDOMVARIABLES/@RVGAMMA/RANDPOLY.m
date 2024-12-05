function p = RANDPOLY(r)

param=getparam(r);
p=POLYLAGUERRE(param.a);
p=setnumber(p,getnumber(r));
