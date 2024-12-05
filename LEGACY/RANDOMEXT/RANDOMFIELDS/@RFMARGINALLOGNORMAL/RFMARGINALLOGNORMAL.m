function ma = RFMARGINALLOGNORMAL(mu,sigma,x0,varargin)
% function ma = RFMARGINALLOGNORMAL(mu,sigma,x0)
% function ma = RFMARGINALLOGNORMAL(mu,sigma,x0,'selfstat')
% 
ma = struct();
ma.mu = mu;
ma.sigma = sigma ; 
ma.x0 = x0 ; 
ma.selfstat = ischarin('selfstat',varargin);

ma = class(ma,'RFMARGINALLOGNORMAL',RFMARGINAL());

