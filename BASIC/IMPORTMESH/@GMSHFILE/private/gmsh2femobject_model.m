function M = gmsh2femobject_model(indim,file,dim,varargin)
% function M = gmsh2femobject_model(indim,file,dim,varargin)

if nargin==2
    dim = indim;
end

fid = fopen(file);
if fid==-1
    error(['Le fichier ' file ' n''existe pas']);
end

[node,elem] = gmsh2femobject_nodeelem(fid,indim);

fclose(fid);

M = cell(1,4);
for i=0:3
    M{i+1} = MODEL(indim);
    M{i+1} = addnode(M{i+1},node);
end

for i=1:length(elem)
    M{getdim(elem{i})+1} = addelem(M{getdim(elem{i})+1},elem{i});
end
for i=1:4
    M{i} = removenodewithoutelem(M{i});
end

if ~ischarin('nofaces',varargin)
%     fprintf('creating facets ridges peaks...\n')
%     if numel(dim)>1
%         error('demander un seul model')
%     end
%     subfun = {@setfacets , @setridges , @setpeaks};
%     for k=dim:-1:1
%         M{dim+1} = subfun{k}(M{dim+1},M{k});
%     end
    
    for i=indim:-1:min(dim)
        M{i+1} = setfacets(M{i+1},M{i});
        if i>1
            M{i+1} = setridges(M{i+1},M{i-1});
            if i>2
                M{i+1} = setpeaks(M{i+1},M{i-2});
            end
        end
    end
    
end

M = M{dim+1};



