function varargout = gmshfile(D,cl,numberpoints,numberlines,numberlineloop,numbersurface,varargin)
% function G = gmshfile(D,cl,numberpoints,numberlines,numberlineloop,numbersurface)
% D : POLYGON
% cl : characteristic length

if nargin<=2
    n = numel(D.P);
    numberpoints = 1:n;
    numberlines = 1:n;
    numberlineloop = n+1;
    numbersurface = 1;
elseif nargin==5
    numbersurface = [];
end

G = GMSHFILE();
P = getvertices(D);
G = createpoints(G,P(1:n),cl,numberpoints);
G = createcontour(G,numberpoints(1:n),numberlines,numberlineloop);
if ~isempty(numbersurface)
    G = createplanesurface(G,numberlineloop,numbersurface);
    if ischarin('recombine',varargin)
        G = recombinesurface(G,numbersurface);
    end
end

varargout{1} = G;
varargout{2} = numbersurface;
