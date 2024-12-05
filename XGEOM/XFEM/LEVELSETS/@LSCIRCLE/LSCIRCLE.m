function ls=LSCIRCLE(cx,cy,r,normtype,varargin)
% function ls=LSCIRCLE(cx,cy,r,normtype,varargin)
% cx, cy : coordonnee du centre
% r : rayon
% normtype : 1 2 ou Inf (2 par defaut)

ls=struct();

param=cell(4,2);
param(:,1)={'cx','cy','r','normtype'};

if nargin<4
    normtype = 2 ;
end

param(:,2)={cx,cy,r,normtype};

lsp = LEVELSET(param,varargin{:});
ls=class(ls,'LSCIRCLE',lsp);
