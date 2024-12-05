function u = solve(K,f)
% function u = solve(K,f)

if isa(K,'FEARRAY')
    rep = K.ddlfree;
    repb = K.ddlbloque;
else
    rep = f.ddlfree;
    repb = f.ddlbloque;
end


u = myzeros(size(f));
if isa(K,'FEARRAY') && isa(f,'FEARRAY')
    u(rep,:) = K.MYDOUBLE(rep,rep)\f.MYDOUBLE(rep,:);
elseif isa(K,'FEARRAY')
    u(rep,:) = K.MYDOUBLE(rep,rep)\f(rep,:);
else
    u(rep,:) = K(rep,rep)\f.MYDOUBLE(rep,:);
end

u = FEVECTOR(u,repb);
