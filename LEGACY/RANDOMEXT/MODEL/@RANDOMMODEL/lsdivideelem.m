function S = lsdivideelem(M)

if length(M.ls)==0
    error('pas de levelset definie pour le modele')
end
if israndom(M.ls)
    error('on ne peut pas diviser les elements pour des levelset aleatoires')
end
M = lseval(M);
M = lssplitelem(M);
S = MODEL(M.mode);
S = addnode(S,M.node);


for i=1:M.nbgroupelem
    elem = M.groupelem{i};
    
    switch getlstype(elem)
        case {'indomain','in','out'}
            S=addelem(S,elem,'norenumelem');
        case 'cut'
            ls = getlevelset(M.ls,getlsnumber(elem));
            [elemcutin,elemcutout,nodeplus]=lsdivideelem(elem,ls,S.node);
            
            S = addnode(S,nodeplus);
            for k=1:length(elemcutin)
                S = addelem(S,elemcutin{k},'norenumelem');
            end
            for k=1:length(elemcutout)
                S = addelem(S,elemcutout{k},'norenumelem');
            end
    end
    
end

M = unique(M);

