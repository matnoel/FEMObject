function RF = RFLOGNORMAL(mu,sigma,x0,C,varargin)
% RF = RFLOGNORMAL(mu,sigma,x0,C,varargin)
% LOGNORMAL STOCHASTIC FIELD 
% mu : mean (can be an inline function)
% sigma : standard deviation (can be an inline function)
% C : correlation 
% mu et sigma, C sont les moyennes, ecart-types, fonction de correlation
% du champ gaussien sous-jacent
%
% loi marginale:
% p(x+x0) = (x*sigma*sqrt(2pi))^-1 * exp(-(log(x)-mu)^2/(2sigma^2))

RF=struct();
type = 'lognormal';

marginal=RFMARGINALLOGNORMAL(mu,sigma,x0,varargin{:});

RF.correl = C;
if ~isa(RF.correl,'RFCORREL')
    error('third argument must be the correlation')
end

RFp = RANDFIELD(marginal,RF.correl);
RF = class(RF,'RFLOGNORMAL',RFp);


