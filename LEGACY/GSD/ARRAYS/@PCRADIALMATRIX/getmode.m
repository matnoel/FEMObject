function u=getmode(u,k)

u.m = length(k);
u.V = MULTIMATRIX(u.V{k},size(u),[u.m,1]);
u.L = getL(u,k);
if ~isempty(u.DLmasse)
u.DLmasse = MULTIMATRIX(u.DLmasse{k},size(u.DLmasse),[u.m,1]);
end
u.D = u.D(k,k);

