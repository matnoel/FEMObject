function ls=LSRECTANGLE(P1x,P1y,P2x,P2y,varargin)
% function ls=LSRECTANGLE(P1x,P1y,P2x,P2y,varargin)
% P1 P2 : 2 opposite corners

ls=struct();

param={'P1x',P1x;'P1y',P1y;'P2x',P2x;'P2y',P2y};

lsp = LEVELSET(param,varargin{:});
ls=class(ls,'LSRECTANGLE',lsp);
