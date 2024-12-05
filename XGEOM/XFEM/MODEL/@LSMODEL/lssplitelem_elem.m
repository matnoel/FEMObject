function [S,I_elemtot,H]=lssplitelem_elem(S,varargin)
% [S,I_elemtot,H]=lssplitelem_elem(S,PCls,varargin)

if isa(S.ls,'LEVELSET')
    S.ls=LEVELSETS(S.ls);
end

H=RANDPOLYS();
for i=1:length(S.ls)
    
    ls = S.ls{i};
    
    if israndom(ls)
        tol = getcharin('tolsplit',varargin,1e-1);
        anisotrope = getcharin('anisotrope',varargin,1);
        statepoint = getcharin('statepoint',varargin,0);
        [S,I_elemtot,H] = lsrandomsplit_elem(ls,tol,anisotrope,statepoint,'MODEL',S,varargin);
    else
        ls = lseval(ls,S);
        S.ls{i}=ls;
        
        for p=1:S.nbgroupelem
            typeenrich = getcharin('lsenrichtype',varargin,1);
            S.groupelem{p}=setlsenrichtype(S.groupelem{p},typeenrich);
            
            
            [elemin,elemcut,elemout] = lssplitelem(S.groupelem{p},ls,S.node);
            %elemin = setmaterial(elemin,S.groupelem{p});
            %elemout = setmaterial(elemout,S.groupelem{p});
            %elemcut = setmaterial(elemin,S.groupelem{p});
            
            enrich = ~ischarin('noenrich',varargin) & ~isempty(getmaterial(ls));
            
            
            elemin = setlsenrich(elemin,0);
            elemin = setlstype(elemin,'in');
            elemin = setlsnumber(elemin,getnumber(ls));
            
            elemcut = setlsenrich(elemcut,enrich);
            elemcut = setlstype(elemcut,'cut');
            elemcut = setlsnumber(elemcut,getnumber(ls));
            
            if isempty(getmaterial(ls))
                elemout = setlstype(elemout,'out');
                %    elemout = setmaterial(elemout,[]);
            end
            
            if isempty(getlsnumber(elemout))
                elemout = setlsnumber(elemout,getnumber(ls));
            end
            
            if ischarin('allenrich',varargin) | (typeenrich>1 & enrich)
                elemin = setlsenrich(elemin,1);
                elemout = setlsenrich(elemout,1);
                elemcut = setlsenrich(elemcut,1);
            end
            
            S.groupelem{p}=elemin;
            S.groupelem{S.nbgroupelem+1}=elemcut;
            S.groupelem{S.nbgroupelem+2}=elemout;
            S.nbgroupelem=length(S.groupelem);
        end
        
        
    end
    S = removeemptygroup(S);
    
    
end


S.H = H;


