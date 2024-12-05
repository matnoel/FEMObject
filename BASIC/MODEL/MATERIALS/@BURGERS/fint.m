function fe = fint(mat,elem,xnode,xgauss,qe,varargin)
% function fe = fint(mat,elem,xnode,xgauss,qe,varargin)

% [se,B] = sigma(mat,elem,xnode,xgauss,qe,varargin{:});

% fe = B'*se;

N = calc_N(elem,xnode,xgauss);
B = calc_B(elem,xnode,xgauss);
u = N*qe;
gradu = B*qe;
k = evalparam(mat,'k',elem,xnode,xgauss);

switch getparam(mat,'formulation')
    case 1
        fe = B'*(k*gradu-1/2*(1-u).*u);
        
    case 2
        alpha = getparam(mat,'alpha');
        Ntilde = N + alpha*B ;
        fe = B'*k*gradu + Ntilde'*(1/2-u)*gradu;
        
end

