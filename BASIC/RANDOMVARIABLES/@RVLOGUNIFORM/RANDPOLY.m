function p = RANDPOLY(r)

param=getparam(r);
p=POLYLOGUNIFORM(r);
p=setnumber(p,getnumber(r));
