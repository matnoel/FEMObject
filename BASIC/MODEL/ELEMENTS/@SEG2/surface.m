function varargout = surface(elem,node,varargin)
% function varargout = surface(elem,node,varargin)

nodecoord = double(getcoord(node));
colnode = getcharin('FaceVertexCData',varargin);
nodecoord = [nodecoord,colnode];
[a,connec] = ismember(getconnec(elem),getnumber(node)) ;

H = patch('faces',connec,'vertices',nodecoord,'Facelighting','phong',varargin{:});

if nargout>=1
    varargout{1} = H;
end
