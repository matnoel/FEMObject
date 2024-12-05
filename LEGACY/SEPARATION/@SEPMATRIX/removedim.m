function uu = removedim(u,dim)
%function uu = removedim(u,dim)

uu.dim=u.dim-numel(dim);
uu.m=u.m;
uu.alpha=u.alpha;

d=1:u.dim;
d(dim)=[];

uu.F=u.F(:,d);

uu=SEPMATRIX(uu);
