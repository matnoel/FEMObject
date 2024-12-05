function varargout = gmshfile(D,cl,numberpoints,numberlines,numberlineloop,numbersurface,varargin)
% function G = gmshfile(D,cl,numberpoints,numberlines,numberlineloop,numbersurface)
% D : QUADRANGLE
% cl : characteristic length

if nargin<=2
    numberpoints = 1:4;
    numberlines = 1:4;
    numberlineloop = 5;
    numbersurface = 1;
elseif nargin==5
    numbersurface = [];
end

G = GMSHFILE();
P = getvertices(D);
G = createpoints(G,P(1:4),cl,numberpoints);
G = createcontour(G,numberpoints(1:4),numberlines,numberlineloop);
if ~isempty(numbersurface)
    G = createplanesurface(G,numberlineloop,numbersurface);
    if ischarin('recombine',varargin)
        G = recombinesurface(G,numbersurface);
    end
end

varargout{1} = G;
varargout{2} = numbersurface;
