function ls=LSSPHERE(cx,cy,cz,r,normtype,varargin)
% function ls=LSSPHERE(cx,cy,cz,r,normtype,varargin)
% cx, cy, cz : coordonnee du centre
% r : rayon
% normtype : 1 2 ou Inf (2 par defaut)

ls=struct();

param=cell(5,2);
param(:,1)={'cx','cy','cz','r','normtype'};

if nargin<5
    normtype = 2 ;
end

param(:,2)={cx,cy,cz,r,normtype};

lsp = LEVELSET(param,varargin{:});
ls=class(ls,'LSSPHERE',lsp);
