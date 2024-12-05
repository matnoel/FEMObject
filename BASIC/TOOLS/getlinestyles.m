function c = getlinestyles(i,varargin)
% function c = getlinestyles(i)
% i : entier
% permet d'obtenir des styles de courbes differents
%  '-','--','-.',':','.'...
%
% function c = getcourbestyles(i,'nomarker')
% uniquement des couleurs differentes


c = {'-','--','-.',':','.'};

i = mod(i-1,length(c))+1;
c = c{i};

