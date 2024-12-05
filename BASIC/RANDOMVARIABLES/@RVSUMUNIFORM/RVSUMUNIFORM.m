function rv = RVSUMUNIFORM(varargin)
% function rv = RVSUMUNIFORM(m,x0,x1)
% variable aleatoire somme de m variables uniformes ind�pendantes 
% bornes [x0,x1]
% 
% function rv = RVSUMUNIFORM(m,mu,sigma,'stat')
% mu et sigma : moyenne et ecart-type de la somme
%
%
% See also RVUNIFORM, RVNORMAL, RVBETA, RVGAMMA, RVLOGNORMAL, RVWEIBULL,
% RVGUMBEL, RVEXTREMEVALUE
rv=struct();
switch nargin
case 0    
    param.m = 1;
    param.x0 = -1;
    param.x1 = 1;
    param.mu = 0;
    param.sigma = 1/sqrt(3);
case 1
    param.m = varargin{1};
    param.x0 = -1/m;
    param.x1 = 1/m;
    param.mu = 0;
    param.sigma = (param.x1-param.x0)/2/sqrt(3*param.m);
case 3
    param.m = varargin{1};
    param.x0 = varargin{2};
    param.x1 = varargin{3};    
    param.mu = (param.x0+param.x1)/2;
    param.sigma = (param.x1-param.x0)/2/sqrt(3*param.m);
case 4
    if ~strcmp(varargin{4},'stat')
        error('argument 4 doit etre ''stat''')
    end
    param.m = varargin{1};
    mu = varargin{2};
    sigma = varargin{3};
    param.x0 = mu-sigma*sqrt(3*param.m);
    param.x1 = mu+sigma*sqrt(3*param.m);    
    param.mu = mu;
    param.sigma = sigma;
otherwise
    error('pas prevu')
end

rvp = RANDVAR('sumunif',param);
rv = class(rv,'RVSUMUNIFORM',rvp);
superiorto('RANDVARS');
