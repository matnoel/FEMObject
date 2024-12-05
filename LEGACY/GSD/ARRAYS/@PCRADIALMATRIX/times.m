function w=times(u,v)

if isa(u,'PCRADIALMATRIX') & ~isa(v,'PCRADIALMATRIX')
    if isa(v,'double')
        w=u;
        w.V = u.V.*v;
        if normest(v)<eps
            w = PCRADIALMATRIX(size(w.V),u.POLYCHAOS);
        end
        
    elseif isa(v,'PCMATRIX')
        if all(size(v)==1)
            w=mtimes(u,v);
        else
            
            u = actualise_ximasse(u);
            w = PCMATRIX([size(u,1),size(v,2)],getPC(v));
            for i=1:u.m
                w = w + multimtimes(u.DLmasse{i},u.V{i}.*v);
            end
        end
        
        %w = times(expand(u),v);
        
    end
    
elseif isa(v,'PCRADIALMATRIX') & ~isa(u,'PCRADIALMATRIX')
    
    if  isa(u,'double')
        w=v;
        w.V = u.*v.V;
        if normest(u)<eps
            w = PCRADIALMATRIX(size(w.V),v.POLYCHAOS);
        end
        
    elseif isa(u,'PCMATRIX')
        
        if all(size(u)==1)
            w=mtimes(u,v);
        else
            %w = times(u,expand(v));
            
            v = actualise_ximasse(v);
            w = PCMATRIX([size(u,1),size(v,2)],getPC(u));
            for i=1:v.m
                w = w + multimtimes(v.DLmasse{i},u.*v.V{i});
            end
        end
        
    end
    
    
    
elseif isa(u,'PCRADIALMATRIX') & isa(v,'PCRADIALMATRIX')
    k=0;
    u=actualise_ximasse(u);
    
    for i=1:u.m
        for j=1:v.m
            k=k+1 ;
            V{k} = u.V{i}.*v.V{j};
            
            L{k} = multimtimes(u.DLmasse{i},v.L(j));
            
        end
    end
    w = PCRADIALMATRIX(V,size(V{1}),L);
    if isdouble(u) & isdouble(v)
        w = cell2mat(w);
    end
    
    
else
    error('times not defined')
end

