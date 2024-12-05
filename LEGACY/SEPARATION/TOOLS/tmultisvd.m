function [u,result] = tmultisvd(b,varargin)
% function [u,result] =  tmultisvd(b,varargin)

if isa(b,'double')
    b=tensor(b);
end

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
errorder = zeros(1,param.maxorder);


switch param.errorindicator
    case 'residual'
        normref = norm(b);
    case 'reference'
        param.reference = b;
        normref = norm(param.reference);
    otherwise
        normref=norm(b);
end


valmin = [];
valmax = [];
errorderinf = zeros(1,param.maxorder);


for i=1:param.maxorder
    if param.display
        fprintf('order #%d \n',i)
    end
    
    U=TSEPMATRIX(tucker_als(bu,ones(1,dim)));
    u=u+U;
    u=orth(u);
    gu=gathervectors(u);
    alpha=ttm(u.alpha,gu.F{1},1);
    for j=2:u.dim
        alpha=ttm(alpha,gu.F{j},j);
    end
    u=splitvectors(gu);
    
    errorder(i)=abs(normref^2-norm(u.alpha)^2)/normref^2;
    
    if param.display
        fprintf('  order #%d - error = %d \n',i,errorder(i))
    end
    
    if errorder(i)<param.tol
        break
    end
    
    bu = b - expand(u);
end

result.error = errorder;
result.errorinf = errorderinf;
result.valmin = valmin ;
result.valmax = valmax;
