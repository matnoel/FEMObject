function mat = BURGERS(varargin)
% function mat = BURGERS('k',k)
% equation de burgers 1/2*(u(1-u)),x - k u,xx = 0

mat = struct();

param = struct(varargin{:});

param.k = getcharin('k',varargin); 
param.formulation = getcharin('formulation',varargin);
param.alpha = getcharin('alpha',varargin,0);

matp = MATERIAL('BURGERS',param);
mat = class(mat,'BURGERS',matp);
