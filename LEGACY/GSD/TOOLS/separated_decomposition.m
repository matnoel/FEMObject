function [u,result] = separated_decomposition(b,varargin)

S = SEPARATEDSOLVER(varargin{:});

[u,result] = separated_decomposition(S,b);

