function varargout = gmshfilewithpoints(C,P,clC,clP,numbercenter,numberpoints,numberlines,numberlineloop,numberembeddedpoints,numbersurface,varargin)
% function G = gmshfilewithpoints(C,P,clC,clP,numbercenter,numberpoints,numberlines,numberlineloop,numberembeddedpoints,numbersurface)
% C : ELLIPSE
% P : POINT
% clC, clP : characteristic length

if ~iscell(P)
    P = {P};
end
if nargin<=3
    clP = clC;
end
if length(clP)==1
    clP = repmat(clP,1,length(P));
end
if nargin<=4
    numbercenter = 5;
    numberpoints = 1:4;
    numberlines = 1:4;
    numberlineloop = 5;
    numberembeddedpoints = 5+(1:length(P));
    numbersurface = 1;
elseif nargin==9
    numbersurface = [];
end

G = GMSHFILE();
P = getvertices(C);
G = createpoint(G,[C.cx,C.cy],clC,numbercenter);
G = createpoints(G,PC,clC,numberpoints);
G = createellipsecontour(G,numbercenter,numberpoints,numberlines,numberlineloop);
if ~isempty(numbersurface)
    G = createplanesurface(G,numberlineloop,numbersurface);
end
G = createpoints(G,P,clP,numberembeddedpoints);
if ~isempty(numbersurface)
    G = embedpointsinsurface(G,numberembeddedpoints,numbersurface);
    if ischarin('recombine',varargin)
        G = recombinesurface(G,numbersurface);
    end
end

varargout{1} = G;
varargout{2} = numbersurface;
