function u = removedim(u,dim)
d=1:u.dim;
d(dim)=[];
u.dim=u.dim-numel(dim);
u.F=u.F(:,d);

