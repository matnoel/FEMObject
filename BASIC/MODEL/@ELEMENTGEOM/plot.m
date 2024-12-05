function varargout = plot(elem,node,varargin)
% function varargout = plot(elem,node,varargin)

if nargin==2
    options = patchoptions(getdim(elem));
elseif getdim(elem)==0
    options = patchoptions(getdim(elem),varargin{:});
else
    options = varargin;
end

nodecoord = double(getcoord(node));
connec = calc_conneclocal(elem,node);

faces = patchfaces(elem,connec);
facevertexcdata = getcharin('facevertexcdata',options);

if ~isempty(facevertexcdata) && size(facevertexcdata,1)==getnbelem(elem)
    facevertexcdata = repmat(facevertexcdata',size(faces,1)/getnbelem(elem),1);
    facevertexcdata = facevertexcdata(:);
    options = setcharin('facevertexcdata',options,facevertexcdata);
end

H = patch('faces',faces,'vertices',nodecoord,options{:});

if nargout>=1
    varargout{1} = H;
end
