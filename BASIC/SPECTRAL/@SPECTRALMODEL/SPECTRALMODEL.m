function x = SPECTRALMODEL(dim,L,varargin)
% function x = SPECTRALMODEL(dim,L)
% Si dim=1 ou 2
% L : POLYLAGRANGE ou POLYFELAGRANGE
if nargin==0
    x.dim = [];
    x.node = NODE();
    x.grid = [];
    x.L = [];
    x.bc = [];   
    x = class(x,'SPECTRALMODEL');
else
    x.dim = dim;
    if dim==1
        x.node = NODE(POINT(getpoints(L)));
        x.grid.X = getpoints(L);
    else
        [X,Y]=meshgrid(getpoints(L),getpoints(L));    
        x.node = NODE(POINT([X(:),Y(:)]));    
        x.grid.X = X;
        x.grid.Y = Y;

    end
    x.L = L;
    x.bc = [];
    x = class(x,'SPECTRALMODEL');
end
