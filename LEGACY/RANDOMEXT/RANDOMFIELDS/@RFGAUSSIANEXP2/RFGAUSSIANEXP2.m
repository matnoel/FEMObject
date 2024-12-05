function RF = RFGAUSSIANEXP2(mu,sigma,L)
% function RF = RFGAUSSIANEXP2(mu,sigma,L)
% GAUSSIAN STOCHASTIC FIELD WITH EXPONENTIAL SQUARE CORRELATION
% mu : mean (can be an inline function)
% sigma : standard deviation (can be an inline function)
% L : correlation length

RFp = RFGAUSSIAN(mu,sigma,EXP2CORREL(L));
RF=struct();
RF = class(RF,'RFGAUSSIANEXP2',RFp);


