function c = getcourbestyles(i,varargin)
% function c = getcourbestyles(i)
% i : entier
% permet d'obtenir des styles de courbes differents
%   'b*-','rs-','ko-','m+-',...
%
% function c = getcourbestyles(i,'nomarker')
% uniquement des couleurs differentes


if ischarin('nomarker',varargin) &&  ischarin('solid',varargin) 
    c = {'k-','b-','r-','m-','g-','c-','y-'};
elseif ischarin('nomarker',varargin) 
    c = {'b-','r--','k-.','m:','g-','c--','y-.'};
else
    c = {'b*-','rs-','ko-','m+-','gd-','cv-','y>-',...
        'b^-','r<-','k>-','mx-','g*-','c<-','yo-',...
        'bs-','ro-','k+-','md-','g^-','c*-','yd-',...
        'b<-','r>-','kx-','m*-','gs-','cx-','ys-','k*-'};
end


i = mod(i-1,length(c))+1;
c = c{i};

