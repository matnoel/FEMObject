function [subgaussin,subgaussout] = calc_lssubgauss(elem,ls,order,varargin)
% function [subgaussin,subgaussout] = calc_lssubgauss(elem,ls,order,varargin)

p = calc_gaussorder(elem,order);

[subgaussin,subgaussout] = elem_lssubgauss(elem,ls,p,varargin{:});

subgaussin = permutegaussND(subgaussin);
subgaussout = permutegaussND(subgaussout);

