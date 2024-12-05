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

mmgoptions = '';
if ischarin('mmgoptions',varargin)
    mmgoptions = [mmgoptions ' ' getcharin('mmgoptions',varargin)];
end
if isempty(strfind(mmgoptions,'-out')) % for compatibility with Matlab version < 9.1 (R2016b)
% if ~contains(mmgoptions,'-out') % for Matlab versions >= 9.1 (R2016b)
    mmgoptions = [mmgoptions ' -out ' getfilemesh(G)];
end

G = appendnodedata(G,q);
% switch dim
%     case 2
%         if isempty(strfind(mmgoptions,'-3dMedit')) % for compatibility with Matlab version < 9.1 (R2016b)
%         % if ~contains(mmgoptions,'-3dMedit') % for Matlab versions >= 9.1 (R2016b)
%             mmgoptions = [mmgoptions ' -3dMedit 1']; % to force mmg2d to produce a .mesh file readable by Gmsh (not anymore compatible with Medit)
%         end
%         G = runfilemmg2d(G,mmgoptions);
%     case 3
%         G = runfilemmg3d(G,mmgoptions);
% end
% G = importmesh(G,dim,varargin{:});

gmshoptions = '';
if ischarin('gmshoptions',varargin)
    gmshoptions = [gmshoptions ' ' getcharin('gmshoptions',varargin)];
end
if isempty(strfind(gmshoptions,'-bmg')) % for compatibility with Matlab version < 9.1 (R2016b)
% if ~contains(gmshoptions,'-bmg') % for Matlab versions >= 9.1 (R2016b)
    gmshoptions = [gmshoptions ' -bmg ' getfilemsh(G) '.pos'];
end
G = remesh(G,dim,'gmshoptions',gmshoptions);

n=max(nargout,1);
varargout = cell(1,n);
[varargout{:}] = gmsh2femobject(indim,G,dim:-1:dim-n+1,varargin{:});
