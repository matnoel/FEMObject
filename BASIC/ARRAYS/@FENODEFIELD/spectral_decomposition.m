function [u,res] = spectral_decomposition(u,varargin)
% function u = spectral_decomposition(u,varargin)
% see also POLYCHAOS/spectral_decomposition

[u.value,res] = spectral_decomposition(u.value,varargin{:});
