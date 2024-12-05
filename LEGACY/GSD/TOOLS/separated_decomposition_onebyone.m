function [u,result] = separated_decomposition_onebyone(b,varargin)

S = SEPARATEDSOLVER(varargin{:});
S = setparam(S,'onebyone',true);

[u,result] = separated_decomposition(S,b);

