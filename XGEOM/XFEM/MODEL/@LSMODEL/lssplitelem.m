function S=lssplitelem(S,varargin)


for i=1:length(S.ls)
    
    ls = S.ls{i};
    
    %if israndom(ls)
    %tol = getcharin('tolsplit',varargin,1e-12);
    %[S,h]=lsrandomsplit(ls,tol,S,varargin{:});
    %H = RANDPOLYS(H,h);
    %else
    
    if ~israndom(ls)
        ls = lseval(ls,S);
        S = setlevelset(S,ls,i);
    end
    
    newgroups = {};
    
    for p=1:getnbgroupelem(S)
        
        if isa(ls,'LSCRACK')
            if ~israndom(ls)
                [elemcut,elembicut,elemin] = lscracksplitelem(getgroupelem(S,p),ls,getnode(S));
            else
                [elemcut,elembicut,elemin] = lsrandomcracksplitelem(getgroupelem(S,p),ls,getnode(S));
            end
            addielem = {elemin,elemcut,elembicut};
            
        else
            
            if ~israndom(ls)
                [elemin,elemcut,elemout] = lssplitelem(getgroupelem(S,p),ls,getnode(S));
            else
                [elemin,elemcut,elemout] = lsrandomsplitelem(getgroupelem(S,p),ls,getnode(S));
            end
            addielem = {elemin,elemcut,elemout};
        end
        
        newgroups = [newgroups , addielem];
        S = setgroupelem(S,p,[]);
    end
    
    for kk=1:length(newgroups)
        S = addgroupelem(S,newgroups{kk});
    end
    %end
    S = removeemptygroup(S);
    
    
end


