function S = addbc(S,varargin)
% function S = addbc(S,varargin)
% rajout de conditions aux limites au modele
% See also BCOND/addbc

S.BCOND = addbc(S.BCOND,varargin{:});

