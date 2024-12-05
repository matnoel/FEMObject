function ls=LSROUGHCIRCLE(cx,cy,vx,vy,r,argr,varargin)
% function ls=LSROUGHCIRCLE(cx,cy,vx,vy,r,argr,varargin)
% D : MODEL (maillage)
% cx, cy : coordonnee du centre du cercle moyen
% vx, vy : vecteur permettant de parametrer le cercle par un angle theta
% theta etant compt� � partir de ce vecteur
% r : function handle donnant le rayon en fonction de l'angle theta
%   r(theta,argr{:}) doit donner la valeur du rayon dans la direction reperee par
%   l'angle theta
% argr : double 1-by-N donnant N arguments supplementaires
% argr peut etre une RANDVARS

ls=struct();

param=cell(6,2);
param(:,1)={'cx','cy','vx','vy','r','argr'};
if nargin<=5
    argr=[];
end
param(:,2)={cx,cy,vx,vy,r,argr};

lsp = LEVELSET(param,varargin{:});
ls=class(ls,'LSROUGHCIRCLE',lsp);
