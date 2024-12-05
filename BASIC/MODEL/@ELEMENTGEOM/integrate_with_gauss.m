function I = integrate_with_gauss(elem,xnode,gauss,fun,varargin)
% function I = integrate_with_gauss(elem,xnode,gauss,fun,varargin)

fun=fcnchk(fun); 
if size(gauss.coord,1)>1
   gauss.coord = permute(gauss.coord,[4,2,3,1]);
   gauss.w     = permute(gauss.w(:),[2,3,4,1]);
end
    
detJ = calc_detJ(elem,xnode,gauss.coord);
funeval = fun(gauss.coord,elem,xnode,varargin{:});

I = sum(gauss.w*abs(detJ)*funeval,4);
