function ke = diff(mat,elem,xnode,xgauss,varargin)
% function ke = diff(mat,elem,xnode,xgauss,varargin)

% if isa(D2,'FENODEFIELD')
% N = calc_N(elem,xnode,xgauss);  
% D2e = localize(elem,D2);
% D2 = N*D2e;
% DD2 = 
% end
D2 = getparam(mat,'D2');
DN = calc_DN(elem,xnode,xgauss);
ke = DN'*D2*DN;
