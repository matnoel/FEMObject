function RF = RFGAUSSIANEXP(mu,sigma,L)
% function RF = RFGAUSSIANEXP(mu,sigma,L)
% GAUSSIAN STOCHASTIC FIELD WITH EXPONENTIAL CORRELATION
% mu : mean (can be an inline function)
% sigma : standard deviation (can be an inline function)
% L : correlation length

RFp = RFGAUSSIAN(mu,sigma,EXPCORREL(L));
RF=struct();
RF = class(RF,'RFGAUSSIANEXP',RFp);
