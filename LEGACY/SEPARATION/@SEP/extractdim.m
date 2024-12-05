function u = extractdim(u,dim)

u.dim=numel(dim);
u.F=u.F(:,dim);

