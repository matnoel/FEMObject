function M = keepgroupelem(M,varargin)
% function M = keepgroupelem(M,groupelem)
% si groupelem est un double -> on garde les groupes d'elements de groupelem
% si groupelem=='out'      : on garde les groupes d'elements avec un lstype='out'
% si groupelem=='in'       : on garde les groupes d'elements avec un lstype='in'
% si groupelem=='cut'      : on garde les groupes d'elements avec un lstype='cut'
% si groupelem=='indomain' : on garde les groupes d'elements avec un lstype='indomain'
% si groupelem=='bicut'    : on garde les groupes d'elements avec un lstype='bicut'
%
% function M = keepgroupelem(M,group1,group2,...)
% on aplique le meme principe a tous les arguments

for i=1:M.nbgroupelem
    garde = 0;
    for j=1:length(varargin)
        if isa(varargin{j},'double') & ismember(i,varargin{j})
            garde=1;
        elseif isa(varargin{j},'char') & strcmp(getlstype(M.groupelem{i}),varargin{j})
            garde=1;
        end
    end
    if garde==0;
        M.groupelem{i}=[];
    end
end

M = removeemptygroup(M);

