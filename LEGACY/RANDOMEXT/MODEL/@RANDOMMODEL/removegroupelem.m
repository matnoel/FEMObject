function M = removegroupelem(M,varargin)
% function M = removegroupelem(M,groupelem)
% si groupelem est un double -> on elimine les groupes d'elements de groupelem
% si groupelem=='out'      : on elimine les groupes d'elements avec un lstype='out'
% si groupelem=='in'       : on elimine les groupes d'elements avec un lstype='in'
% si groupelem=='cut'      : on elimine les groupes d'elements avec un lstype='cut'
% si groupelem=='indomain' : on elimine les groupes d'elements avec un lstype='indomain'
%
% function M = removegroupelem(M,group1,group2,...)
% on aplique le meme principe a tous les arguments

for j=1:length(varargin)
    if isa(varargin{j},'double')
        M.groupelem(varargin{j})={[]};
    elseif isa(varargin{j},'char')
        
        for i=1:M.nbgroupelem
            if strcmp(getlstype(M.groupelem{i}),varargin{j})
                M.groupelem{i}=[];
            end
        end
        
    end
end

M = removeemptygroup(M);