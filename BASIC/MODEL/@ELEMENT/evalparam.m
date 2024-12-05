function varargout = evalparam(elem,mat,paramname,xnode,xgauss)
% function varargout = evalparam(elem,mat,paramname,xnode,xgauss)

param = getparam(mat,paramname);

if isa(param,'FENODEFIELD')
    [val,x] = evalparam(mat,paramname,elem,xnode,xgauss);
% elseif isa(param,'double') || isa(param,'MYDOUBLE')
%     x = calc_x(elem,xnode,xgauss);
%     val = evalparam(mat,paramname,x);
% elseif isa(param,'inline') || isa(param,'function_handle')
%     x = calc_x(elem,xnode,xgauss);
%     val = evalparam(mat,paramname,x);
else
    x = calc_x(elem,xnode,xgauss);
    val = evalparam(mat,paramname,x);
end

varargout{1} = val;
if nargout>1
    varargout{2} = x ;
end

