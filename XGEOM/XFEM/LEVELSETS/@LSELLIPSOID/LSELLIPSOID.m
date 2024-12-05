function ls=LSELLIPSOID(cx,cy,cz,a,b,c,vax,vay,vaz,vbx,vby,vbz,varargin)
% function ls=LSELLIPSOID(cx,cy,cz,a,b,c,vax,vay,vaz,vbx,vby,vbz,varargin)
% D : MODEL (maillage)
% cx, cy , cz : coordonnee du centre
% a , b , c: longueur des demi-axes
% vax, vay ,vaz : vecteur indiquant la direction de l'axe a
% vbx, vby ,vbz : vecteur indiquant la direction de l'axe b

ls=struct();

param=cell(12,2);
param(:,1)={'cx','cy','cz','a','b','c','vax','vay','vaz','vbx','vby','vbz'};

param(:,2)={cx,cy,cz,a,b,c,vax,vay,vaz,vbx,vby,vbz};

lsp = LEVELSET(param,varargin{:});
ls=class(ls,'LSELLIPSOID',lsp);
