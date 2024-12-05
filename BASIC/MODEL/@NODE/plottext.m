function plottext(u,text,varargin)
% function plottext(u,text)
% text : vecteur de caracteres ou de double
%    contenant une info a afficher pour chaque noeud
%
% function plottext(u,text,'color',color,'fontsize',fontsize)

color = getcharin('color',varargin,'b');
fontsize = getcharin('fontsize',varargin,14); 
plottext(u.POINT,text,'fontsize',fontsize,'color',color);

