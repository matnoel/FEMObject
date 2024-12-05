function S = setlevelsets(S,ls,varargin)
% function S = setlevelsets(S,ls,varargin)

if nargin==2
    S.ls = LEVELSETS(ls);
else
    error('pas prevu')
end
