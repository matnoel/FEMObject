function c = getpointstyles(i,varargin)
% function c = getpointstyles(i)
% i : entier
% permet d'obtenir des styles de courbes differents
%   'b*','rs','ko','m+',...
%

c = {'b*','rs','ko','m+','gd',...
    'b^','r<','k>','mx','g*',...
    'bs','ro','k+','md','g^',...
    'b<','r>','kx','m*','gs'};

i = mod(i-1,length(c))+1;
c = c{i};

