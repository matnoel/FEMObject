function ls=LSROUGHLINE(cx,cy,vx,vy,a,arga,varargin)
% function ls=LSROUGHLINE(cx,cy,vx,vy,a,arga,varargin)
% D : MODEL (maillage)
% cx, cy : coordonnee d'un point de la ligne moyenne
% vx, vy : vecteur indiquant le vecteur normal � la ligne moyenne (normale sortante)
% a : function handle telle que a(xi,arga) donne la distance normale 
%     (signee) entre la surface et la ligne moyenne 
%     xi etant l'abscisse curviligne compt�e � partir du point (cx,cy)
% arga : double 1-by-N donnant des parametres supplementaires a la fonction a
% arga peut etre une RANDVARS

ls=struct();

param=cell(6,2);
param(:,1)={'cx','cy','vx','vy','a','arga'};
if nargin<=5
    arga=[];
end
param(:,2)={cx,cy,vx,vy,a,arga};

lsp = LEVELSET(param,varargin{:});
ls=class(ls,'LSROUGHLINE',lsp);
