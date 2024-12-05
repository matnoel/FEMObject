function numelem=getnumelem(M,groups,choix)

if nargin==1
    groups = 1:M.nbgroupelem;
    numelem = getnumelem(M,groups);
elseif nargin==2
    numelem=[];
    for i=1:length(groups)
        num=getnumber(M.groupelem{groups(i)});
        numelem=[numelem;num(:)];
    end
elseif nargin==3 && strcmp(choix,'local')
    numelem = getnumelem(M,groups);
    [temp,numelem] = ismember(numelem,getnumelem(M));
end