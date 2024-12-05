function I = integrate(elem,xnode,order,fun,varargin)
% function I = integrate(elem,xnode,order,fun,varargin)

fun=fcnchk(fun); 

if isa(xnode,'NODE')
    xnode = xnode(elem);
end

gauss = calc_gauss(elem,order);

detJ = calc_detJ(elem,xnode,gauss.coord);
funeval = fun(gauss.coord,elem,xnode,varargin{:});

I = sum(gauss.w*abs(detJ)*funeval,4);
