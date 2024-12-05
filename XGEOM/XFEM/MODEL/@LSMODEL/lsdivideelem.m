function S = lsdivideelem(M)
% function S = lsdivideelem(M)

if length(M.ls)==0
    error('pas de levelset definie pour le modele')
end
if israndom(M.ls)
    error('on ne peut pas diviser les elements pour des levelset aleatoires')
end
M = lseval(M);
M = lssplitelem(M);
S = MODEL(getmode(M));
S = addnode(S,getnode(M));

for i=1:getnbgroupelem(M)
    elem = getgroupelem(M,i);
    
    switch getlstype(elem)
        case {'indomain','in','out'}
            S = addelem(S,elem,'norenumelem');
            
        otherwise
            ls = getlevelset(M.ls,getlsnumber(elem));
            if isa(ls,'LSCRACK')
                [elemcutin,elemcutout,nodeplus] = lscrackdivideelem(elem,ls,S.node);
            else
                [elemcutin,elemcutout,nodeplus] = lsdivideelem(elem,ls,S.node);
            end
            
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
