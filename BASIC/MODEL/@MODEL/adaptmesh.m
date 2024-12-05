function varargout = adaptmesh(M,q,filename,varargin)
% function varargout = adaptmesh(M,q,filename,varargin)

if ~isa(M,'MODEL')
    varargout = cell(1,nargout);
    [varargout{:}] = adaptmesh(q,M,filename,varargin{:});
    return
end

if israndom(q)
    error('la solution est aleatoire')
end
q = unfreevector(M,q);

indim = getindim(M);
dim = getdim(M);

G = GMSHFILE();
if nargin>=3 && ischar(filename)
    G = setfile(G,filename);
end

% G = exportmesh(G,dim,varargin{:});
% if indim<3
%     file = getfilemesh(G);
%     fid = fopen(file,'r+');
%     while ~strcmp(fgetl(fid),' Dimension')
%     end
%     fseek(fid,0,'cof');
%     fprintf(fid,' %u',indim);
%     fclose(fid);
% end
% G = writefilesol(G,dim,q);
% options = [' -sol ' getfilesol(G)];

options = '';
if ischarin('mmgoptions',varargin)
    options = [options ' ' getcharin('mmgoptions',varargin)];
end
if isempty(strfind(options,'-out')) % for compatibility with Matlab version < 9.1 (R2016b)
% if ~contains(options,'-out') % for Matlab versions >= 9.1 (R2016b)
    options = [options ' -out ' getfilemesh(G)];
end

G = appendnodedata(G,q);
switch dim
    case 2
        if isempty(strfind(options,'-3dMedit')) % for compatibility with Matlab version < 9.1 (R2016b)
            % if ~contains(options,'-3dMedit') % for Matlab versions >= 9.1 (R2016b)
            options = [options ' -3dMedit 1']; % to force mmg2d to produce a .mesh file readable by Gmsh (not anymore compatible with Medit)
        end
        G = runfilemmg2d(G,options);
    case 3
        G = runfilemmg3d(G,options);
end
G = importmesh(G,dim,varargin{:});

n=max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,dim:-1:dim-n+1,varargin{:});
