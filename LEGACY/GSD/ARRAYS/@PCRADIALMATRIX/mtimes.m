function w=mtimes(u,v)
if isa(v,'logical')
    v=double(v);
end

if isa(u,'logical')
    u=double(u);
end



if isa(u,'PCRADIALMATRIX') && ~isa(v,'PCRADIALMATRIX')
    if isa(v,'PCMATRIX')
        w=u;
        if all(size(v)==1)
            %w.L = mtimes(w.L,v);
            %PCu=getPC(u);
            %for i=1:u.m
            %    L{i} = PCMATRIX(multimtimesold(v,u.DLmasse{i}),[1,1],PCu);
            %end
            %w.L = vertcat(L{:});
            u=actualise_ximasse(u);
            
            w.L = multimtimes(switchmulti(u.DLmasse),v);
            w.L = PCMATRIX(w.L,getPC(v));
            w.POLYCHAOS = getPC(w.L);
        else
            
            u = actualise_ximasse(u);
            if all(size(u)==1)
                ss = size(v);
            else
                ss=[size(u,1),size(v,2)];
            end
            w = PCMATRIX(ss,getPC(v));
            
            for i=1:u.m
                temp = multimtimes(u.DLmasse{i},getmultimatrix(u.V{i}*v));
                w = w + PCMATRIX(temp,ss,getPC(v));
            end
            
            %temp = (double(v)*u.DLmasse);
            %temp = reshape(temp,size(v,1),size(v,2)*size(temp,2));
            
            %temp = u.V * temp;
            %w=multisum(temp);
            %w=PCMATRIX(w,[size(u,1),size(v,2)],getPC(u));
            
        end
        
    elseif isa(v,'double')
        
        w=u;
        w.V = u.V*v;
        try
            nn = normest(v(:));
        catch
            nn = norm(v);
        end
        if nn<eps
            w = PCRADIALMATRIX(size(w.V),u.POLYCHAOS);
        end
        
    end
    
elseif isa(v,'PCRADIALMATRIX') && ~isa(u,'PCRADIALMATRIX')
    
    if isa(u,'PCMATRIX')
        w=v;
        if all(size(u)==1)
            %w.L = mtimes(u,w.L);
            %PCv=getPC(v);
            %for i=1:v.m
            %    L{i} = PCMATRIX(multimtimesold(u,v.DLmasse{i}),[1,1],PCv);
            %end
            %w.L = vertcat(L{:});
            
            v=actualise_ximasse(v);
            w.L = multimtimes(switchmulti(v.DLmasse),u);
            w.L = PCMATRIX(w.L,getPC(u));
            w.POLYCHAOS = getPC(w.L);
            
        else
            %w=mtimes(u,expand(v));
            %v = calc_ximasse(v);
            %temp = v.DLmasse * double(u).' ;
            %temp = reshape(temp,size(temp,1)*size(u,1),size(u,2));
            %temp = temp*v.V;
            %w=multisum(temp).';
            %w=PCMATRIX(w,[size(u,1),size(v,2)],getPC(u));
            %w = u*v.V;
            v = actualise_ximasse(v);
            if all(size(v)==1)
                ss = size(u);
            else
                ss=[size(u,1),size(v,2)];
            end
            w = PCMATRIX(ss,getPC(u));
            for i=1:v.m
                temp =multimtimes(v.DLmasse{i},getmultimatrix(u*v.V{i}));
                w = w + PCMATRIX(temp,ss,getPC(u));
            end
            
            %norm(w-w2)/norm(w2)
            %norm(w2-w3)/norm(w2)
            
        end
        
    elseif isa(u,'double')
        
        w=v;
        w.V = u*v.V;
        
        try
            nn = normest(u(:));
        catch
            nn = norm(u);
        end
        if nn<eps
            w = PCRADIALMATRIX(size(w.V),v.POLYCHAOS);
        end
        
    end
    
elseif isa(u,'PCRADIALMATRIX') && isa(v,'PCRADIALMATRIX')
    %if all(size(u)==1)
    
    
    k=0;
    u = actualise_ximasse(u);
    V = cell(1,u.m*v.m);
    for i=1:u.m
        for j=1:v.m
            k=k+1 ;
            
            V{k} = u.V{i}*v.V{j};
            
            %L{k} = multimtimesold(v.L(j),u.DLmasse{i});
            L{k} = multimtimes(u.DLmasse{i},v.L(j));
            
        end
    end
    
    w = PCRADIALMATRIX(V,size(V{1}),L);
    if isdouble(u) && isdouble(v)
        w = cell2mat(w);
    end
    
else
    error('mtimes not defined')
end

