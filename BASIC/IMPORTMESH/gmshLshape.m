function varargout = gmshLshape(cl,filename,indim,varargin)
% function varargout = gmshLshape(cl,filename,indim)
% cl : characteristic length
% filename : file name (optional)
% indim : space dimension (optional, 2 by default)

if nargin<3 || isempty(indim)
    indim = 2;
end

P{1} = [0,0,0];
P{2} = [1,0,0];
P{3} = [1,1,0];
P{4} = [2,1,0];
P{5} = [2,2,0];
P{6} = [0,2,0];

G = GMSHFILE();
if nargin>=2 && ischar(filename)
    G = setfile(G,filename);
end

G = createpoints(G,P,cl,1:6);
G = createcontour(G,1:6,1:6,7);
G = createplanesurface(G,7,1);
if ischarin('recombine',varargin)
    G = recombinesurface(G,1);
end
G = createphysicalsurface(G,1,1);
varargin = delonlycharin('recombine',varargin);

n=max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,2:-1:2-n+1,varargin{:});
