function varargout = surface(elem,node,varargin)

if nargin==2
    options = patchoptions(getdim(elem));
elseif getdim(elem)==0
    options = patchoptions(getdim(elem),varargin{:});
else
    options = varargin;
end

nodecoord = double(getcoord(node));
connec = calc_conneclocal(elem,node) ;

faces = patchfaces(elem,connec);
facevertexcdata = getcharin('facevertexcdata',options);

if ~isempty(facevertexcdata) && size(facevertexcdata,1)==getnbelem(elem)
    facevertexcdata = repmat(facevertexcdata',size(faces,2),1);
    facevertexcdata = facevertexcdata(:);
    options = setcharin('facevertexcdata',options,facevertexcdata);
    con = faces';
    nodecoord = nodecoord(con(:),:);
    nodecoord = [nodecoord,facevertexcdata(:)];
    con(:) = 1:numel(faces);
    faces = con';
    
elseif ~isempty(facevertexcdata) && size(facevertexcdata,1)==size(nodecoord,1)
    nodecoord = [nodecoord,facevertexcdata];
end
options = delonlycharin('surface',options);

H = patch('faces',faces,'vertices',nodecoord,options{:});

if nargout>=1
    varargout{1} = H;
end
