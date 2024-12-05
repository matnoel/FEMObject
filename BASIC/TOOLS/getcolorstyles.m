function c = getcolorstyles(i,varargin)
% function c = getcolorstyles(i)
% i : entier
% permet d'obtenir des styles de courbes differents
%   'b','r','k','m',...
%

c = {'b','r','y','g','m','c','k'};

i = mod(i-1,length(c))+1;
c = c{i};


