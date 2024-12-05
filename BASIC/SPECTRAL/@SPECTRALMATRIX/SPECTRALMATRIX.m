function x = SPECTRALMATRIX(a,L,varargin)
% function x = SPECTRALMATRIX(a,L)
% L : POLYLAGRANGE ou POLYFELAGRANGE
% a : double de taille getnbpoints(L)

if nargin==0 
    x.value = [];
    x.L = L;
    x = class(x,'SPECTRALMATRIX');
end
