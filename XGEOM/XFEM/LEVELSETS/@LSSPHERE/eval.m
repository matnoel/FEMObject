function lsx=eval(ls,x,cx,cy,cz,r,normtype)

cx= MYDOUBLEND(reshape(full(cx),[1,1,1,numel(cx)]));
cy= MYDOUBLEND(reshape(full(cy),[1,1,1,numel(cy)]));
cz= MYDOUBLEND(reshape(full(cz),[1,1,1,numel(cz)]));
r = MYDOUBLEND(reshape(full(r),[1,1,1,numel(r)]));

lsx=distance(POINT(x),POINT([cx,cy,cz]),normtype)-r;
