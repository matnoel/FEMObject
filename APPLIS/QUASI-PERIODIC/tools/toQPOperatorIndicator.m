function x = toQPOperatorIndicator(x,strict)
% x = toQPOperatorIndicator(x,strict)
% Variant of TuckerLikeTensor/toOperator for QUASI-PERIODIC library.
% Acts as toOperator on the first x.order-1 dimensions. Along the last
% dimension, every vector v is transformed into v*ones(1,size(v,1)); if 
% boolean 'strict' is true, then it is transformed into v*v'.
% The goal is to reproduce the matrix operation
% M(~x,:) = 0; 
% and, if 'strict' is true,
% M(:,~x) = 0;
% where M is a matrix and x is a boolean-like column (an `indicator'), by 
% X.*MT;
% where MT is the (QP) tensor representation of M and X =
% toQPOperatorIndicator(x).

% Mesoscopic dimensions (same code as toOperator)
x.space = TSpaceOperators(x.space);
x.space.spaces(1:x.order-1) = cellfun(@(y) ...
    cellfun(@(z)spdiags(z,0,size(z,1),size(z,1)),y,'uniformoutput',false),...
    x.space.spaces(1:x.order-1),'uniformoutput',false);

% Microscopic dimension
if strict
    fun = @(z) z*z';
else
    fun = @(z) z*ones(1,size(z,1));
end
x.space.spaces{x.order} = cellfun(fun,x.space.spaces{x.order},...
    'UniformOutput',false) ;
x.space =  updateProperties(x.space);
x =  updateProperties(x);
end