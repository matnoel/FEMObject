function RF = RFGAUSSIAN(varargin)
% function RF = RFGAUSSIAN(mu,sigma,C)
% GAUSSIAN STOCHASTIC FIELD WITH CORRELATION C
% mu : mean (can be an inline function)
% sigma : standard deviation (can be an inline function)
% C : correlation 
RF=struct();
type = 'gaussian';
RF.mu=varargin{1};
RF.sigma = varargin{2};
RF.correl = varargin{3};
marginal=RFMARGINAL(RVNORMAL,'mu',RF.mu,'sigma',RF.sigma);
if ~isa(RF.correl,'RFCORREL')
    error('third argument must be the correlation')
end

RFp = RANDFIELD(marginal,RF.correl);
RF = class(RF,'RFGAUSSIAN',RFp);
