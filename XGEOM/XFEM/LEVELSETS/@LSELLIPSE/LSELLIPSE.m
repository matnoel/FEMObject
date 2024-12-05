function ls=LSELLIPSE(cx,cy,a,b,vx,vy,varargin)
% function ls=LSELLIPSE(cx,cy,a,b,vx,vy,varargin)
% D : MODEL (maillage)
% cx, cy : coordonnee du centre
% a , b : longueur des demi-axes
% vx, vy : vecteur indiquant la direction de l'axe a

ls=struct();

param=cell(6,2);
param(:,1)={'cx','cy','a','b','vx','vy'};

param(:,2)={cx,cy,a,b,vx,vy};

lsp = LEVELSET(param,varargin{:});
ls=class(ls,'LSELLIPSE',lsp);
