function M = concatgroupelem(M,numgroups,varargin)
%function M = concatgroupelem(M)
% concatenation des groupes d'elements identiques
%
%function M = concatgroupelem(M,numgroups)
% concatenation des groupes numgroups
%
%function M = concatgroupelem(M,numgroups,'enforce')
% concatenation des groupes numgroups 
% sans regarder si c'est le meme materiau ou le meme etat levelset

if nargin==1 || isempty(numgroups)
    numgroups=1:M.nbgroupelem;
end
nonumgroups = setdiff(1:M.nbgroupelem,numgroups);
elimgroups = [];
for i=numgroups
 
    mati = getmaterial(M.groupelem{i});
for j=setdiff(1:i-1,nonumgroups)
    matj = getmaterial(M.groupelem{j});
    if isempty(mati) && isempty(matj)
        okmat=1;
    else
        okmat = (mati == matj);
    end
    okmat = okmat | ischarin('enforce',varargin);
    okls = lsdatacmp(M.groupelem{j},M.groupelem{i}) | ischarin('enforce',varargin) ;
    
if strcmp(class(M.groupelem{i}),class(M.groupelem{j})) && okmat && okls
    M.groupelem{j}=concat(M.groupelem{j},M.groupelem{i});
    elimgroups = union(elimgroups,i);
end
    
end
end

M.groupelem(elimgroups) = [];
M.nbgroupelem = length(M.groupelem);
M=changeelemnumber(M);

