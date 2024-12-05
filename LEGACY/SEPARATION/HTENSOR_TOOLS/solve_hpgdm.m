function [u, result] = solve_hpgdm(A,b,HPGD,varargin)

if ischarin('tol',varargin)
    tol = getcharin('tol',varargin);
else
    tol = 1e-8;
end

if ischarin('maxiter',varargin)
    maxiter = getcharin('maxiter',varargin);
else
    maxiter = 20;
end


bu = b;
optsbu.tol = 1e-10;
normref = norm(b);
errorder = zeros(1,maxiter);
v = cell(1,maxiter);

for i=1:maxiter
    v{i} = solve_hpgd_prog(A,bu,HPGD);
    Av = apply_mat_to_vec(A,v{i});
    optsbu.max_rank = max(rank(Av));
    Av = truncate(Av,optsbu);
    
    bu = bu - Av;
    optsbu.max_rank = max(rank(bu));
    bu = truncate(bu,optsbu);
    
    if i==1
        u = v{i};
    else
        u = u+v{i};
    end
    errorder(i) = norm(bu)/normref;
    
    fprintf(' HPGDm - order #%d - error = %d\n',i,errorder(i));
    if errorder(i) < tol
        break
    end
end

result.error = errorder;
result.v = v;


end