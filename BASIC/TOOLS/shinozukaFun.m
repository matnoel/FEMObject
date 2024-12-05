function [V,order] = shinozukaFun(x,Z,Phi,k,c,varargin)
% function [V,order] = shinozukaFun(x,Z,Phi,k,c,varargin)
% Computes a realization of nU Gaussian random fields V evaluated at
% the points in the nx-by-d matrix x with a given spatial correlation structure 
% using the spectral representation of Shinozuka with 
% the 1-by-order-by-nU arrays of values Z and Phi, and 
% the d-by-order^(1/d) arrays of values k and c.
% Inputs:
% x is a d-column matrix that specifies the points where the Gaussian 
% random fields are to be evaluated.
% Z is a 1-by-order-by-nU array of double specifying a realization of 
% random amplitude Z with values in [0,+Inf[.
% Phi is a 1-by-order-by-nU array of double specifying a realization of 
% random phase shift Phi uniformly distributed on [0,2*pi].
% k is a d-by-order^(1/d) array of double specifying the discretization of 
% wave number k in [-pi/lcorr(j),pi/lcorr(j)] for each spatial dimension j.
% c is a d-by-order^(1/d) array of double specifying the discretization of 
% wave number k in [-pi/lcorr(j),pi/lcorr(j)] for each spatial dimension j.
% Outputs:
% V is the nx-by-nU array of double containing a realization of 
% the nU Gaussian random fields evaluatued at the nx points.
% order is the d-dimensional order (number of terms) of the spectral representation.
%
% V = shinozukaFun(x,Z,Phi,k,c) returns the nx-by-nU array of values V
% where the realizations of the nU Gaussian random fields are to be evaluated.
%
% [V,order] = shinozukaFun(...) also returns the order (number of terms) 
% of the spectral representation.

nx = size(x,1);
order = size(Z,2);
nU = size(Z,3);
V = zeros(nx,nU);
for i=1:nU
    z = Z(:,:,i);
    phi = Phi(:,:,i);
    V(:,i) = shinozuka_cpp(c,z,phi,k,x);
end

end
