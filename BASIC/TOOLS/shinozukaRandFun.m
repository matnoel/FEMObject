function [V,order] = shinozukaRandFun(x,Z,Phi,K,C,varargin)
% function [V,order] = shinozukaRandFun(x,Z,Phi,K,C,varargin)
% Computes a realization of nU Gaussian random fields V evaluated at
% the points in the nx-by-d matrix x with a given spatial correlation structure 
% using the randomized spectral representation of Shinozuka with 
% the 1-by-order-by-nU arrays of values Z and Phi,  
% the d-by-order-by-nU array of values K, and 
% the 1-by-order-by-nU array of values C.
% Inputs:
% x is a d-column matrix that specifies the points where the Gaussian 
% random fields are to be evaluated.
% Z is a 1-by-order-by-nU array of double specifying a realization of 
% random amplitude Z with values in [0,+Inf[.
% Phi is a 1-by-order-by-nU array of double specifying a realization of 
% random phase shift Phi uniformly distributed on [0,2*pi].
% K is a d-by-order-by-nU array of double specifying a realization of 
% random wave number K uniformly distributed on [-pi/lcorr(j),pi/lcorr(j)] for each spatial dimension j.
% C is a 1-by-order-by-nU array of double.
% Outputs:
% V is the nx-by-nU array of double containing a realization of 
% the nU Gaussian random fields evaluatued at the nx points.
% order is the order (number of terms) of the spectral representation.
%
% V = shinozukaRandFun(x,Z,Phi,K,C) returns the nx-by-nU array of values V
% where the realizations of the nU Gaussian random fields are to be evaluated.
%
% [V,order] = shinozukaRandFun(...) also returns the order 
% (number of terms) of the spectral representation

nx = size(x,1);
order = size(Z,2);
nU = size(Z,3);
V = zeros(nx,nU);
for i=1:nU
    z = Z(:,:,i);
    phi = Phi(:,:,i);
    k = K(:,:,i);
    c = C(:,:,i);
    V(:,i) = shinozukaRand_cpp(c,z,phi,k,x);
end

end
