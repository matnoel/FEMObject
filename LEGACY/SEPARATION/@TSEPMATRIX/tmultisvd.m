function [u,result] = tmultisvd(b,varargin)
% function [u,result] =  tmultisvd(b,varargin)

dim =length(size(b));
n = size(b);

if ischarin('SEPSOLVER',varargin)
    solver = getcharin('SEPSOLVER',varargin);
else
    solver = SEPSOLVER(dim,varargin{:});
end
param = getparam(solver);

bu = b;
u = TSEPMATRIX(dim);

erriter = zeros(1,param.maxiter);
if ~isfield(param,'maxorder')
    param.maxorder = max(n);
end
errorder = ones(1,param.maxorder);


switch param.errorindicator
    case 'residual'
        normref = norm(b.alpha);
    case 'reference'
        if isempty(param.reference)
            param.reference = b;
        end
        if isa(param.reference,'double')
            normref = norm(param.reference(:));
        else
            normref=norm(param.reference);
        end
    otherwise
        normref=norm(b.alpha);
end


valmin = [];
valmax = [];
errorderinf = zeros(1,param.maxorder);


for i=1:param.maxorder
    if param.display
        fprintf('order #%d \n',i)
    end
    
    U=tucker_als(bu,ones(1,dim));
    
    u=orth(u+U);
    
    gu=gathervectors(u);
    alpha=ttm(u.alpha,gu.F{1},1);
    for j=2:u.dim
        alpha=ttm(alpha,gu.F{j},j);
    end
    u=splitvectors(gu);
    
    bu = orth(b - u);
    
    
    % errorder(i)=abs(normref^2-norm(u.alpha)^2)/normref^2;
    errorder(i)=norm(bu)/normref;
    
    if param.display
        fprintf('  order #%d - error = %d \n',i,errorder(i))
    end
    
    if errorder(i)<param.tol
        break
    end
end

result.error = errorder;
result.errorinf = errorderinf;
result.valmin = valmin ;
result.valmax = valmax;

