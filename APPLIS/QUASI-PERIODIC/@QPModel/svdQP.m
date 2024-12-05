function y = svdQP(model,m)
% y = svdQP(model,m)
% Overload of tool function svdQP, as a QPModel method.
% See tools/svdQP documentation.

y = svdQP(m,getCellNum(model),getTolSVD(model)) ;

end

