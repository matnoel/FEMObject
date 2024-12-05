function M = actualisematerials(M,mat)
% function M = actualisematerials(M)
% numerotation des materiaux non numerotes
%
% function M = actualisematerials(M,mat)
% actualisation des materiaux des elements
% mat : MATERIALS (nouveaux materiaux)

if nargin==1
    for k=1:M.nbgroupelem
        mat{k} = getmaterial(getgroupelem(M,k));
        if isempty(mat{k})
            error(['Pas de materiau defini pour le groupe d''element ' num2str(k)])
        end
    end
    for k=1:length(M.ls)
        mat{k+M.nbgroupelem} = getmaterial(M.ls{k});
    end
    
    for k=1:M.nbgroupelem
        if isempty(getnumber(mat{k}))
            mat{k} = setnumber(mat{k},newnumber(MATERIALS(mat{:})));
            M = setfieldgroupelem(M,'material',mat{k},k);
        end
    end
    for k=1:length(M.ls)
        if isa(mat{k+M.nbgroupelem},'MATERIAL') && isempty(getnumber(mat{k+M.nbgroupelem}))
            mat{k+M.nbgroupelem} = setnumber(mat{k+M.nbgroupelem},newnumber(MATERIALS(mat{:})));
            M.ls{k} = setmaterial(M.ls{k},mat{k+M.nbgroupelem});
        end
    end
    
else
    if ~isa(mat,'MATERIALS')
        error('rentrer un MATERIALS')
    end
    
    for k=1:M.nbgroupelem
        n = getmaterialnumber(getgroupelem(M,k));
        if isempty(n)
            error('les materiaux des elements doivent etre numerotes')
        end
        
        if ismember(n,mat)
            M = setfieldgroupelem(M,'material',getmaterial(mat,n),k);
        end
    end
    
    for k=1:length(M.ls)
        matls = getmaterial(M.ls{k});
        if isa(matls,'MATERIAL')
            n = getnumber(matls);
            if isempty(n)
                error('les materiaux des elements doivent etre numerotes')
            end
            
            if ismember(n,mat)
                M.ls{k} = setmaterial(M.ls{k},getmaterial(mat,n));
            end
        end
    end
    
end

end

