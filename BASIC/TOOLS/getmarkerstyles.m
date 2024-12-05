function c = getmarkerstyles(i,varargin)
% function c = getmarkerstyles(i)
% i : entier
% permet d'obtenir des styles de courbes differents
%   '*','s','o','+',...
%
% function c = getcourbestyles(i,'nomarker')
% uniquement des couleurs differentes


c = {'*','s','o','+','>','<','^'};

i = mod(i-1,length(c))+1;
c = c{i};

