function v = unfreevector(BC,u,k)
% function v = unfreevector(BC,u,k)

if getnbddl(BC)==0
    v=u;
    return
end
if nargin<3
    k=1;
end
if size(u,k)==BC.nbddl
    v=u;
    return
end

if all(gettypes(BC)==0) && isa(u,'double')
    
    rep = cell(ndims(u),1);
    rep(:)={':'};
    rep(k)={BC.ddlfree};
    s=size(u);
    s(k)=BC.nbddl;
    
    if issparse(u)
        v=sparse(s(1),s(2));
    else
        v=zeros(s);
    end
    
    v(rep{:})=u;
    
    if ~all(ishomogeneous(BC)) && all(gettypes(BC)==0)
        
        if k==1
            for i=1:length(BC.BC)
                if s(2)==size(BC.BC{i}.value,2)
                    v(BC.BC{i}.ddlbloque,:)=BC.BC{i}.value;
                else
                    v(BC.BC{i}.ddlbloque,:)=repmat(BC.BC{i}.value,1,s(2));
                end
            end
        elseif k==2
            for i=1:length(BC.BC)
                if s(1)==size(BC.BC{i}.value,1)
                    v(:,BC.BC{i}.ddlbloque)=BC.BC{i}.value;
                else
                    v(:,BC.BC{i}.ddlbloque)=repmat(BC.BC{i}.value,s(1),1);
                end
            end
        else
            error('')
        end
        
    end
    
    
    
elseif all((gettypes(BC)==1) | (gettypes(BC)==0) ) && isa(u,'double')
    s=size(u);
    s(k)=BC.nbddl;
    rep = cell(ndims(u),1);
    rep(:)={':'};
    rep(k)={BC.ddlfree};
    
    if issparse(u)
        v=sparse(s(1),s(2));
    else
        v=zeros(s);
    end
    
    v(rep{:})=u;
    d1=rep;
    d2=rep;
    for i=1:length(BC.BC)
        if BC.BC{i}.type==1
            d1{k}=BC.BC{i}.ddlfree;
            d2{k}=BC.BC{i}.ddlbloque;
            % % START DEBUG % %
            [d1bloque,id1,~] = intersect(d1{k},BC.ddlbloque);
            if ~isempty(d1bloque) % these dof are blocked by another BC
                for j = setdiff(1:length(BC.BC),i)
                    % Find whether BC #j block some of these, and their indices
                    [match,idfree] = ismember(d1bloque,BC.BC{j}.ddlbloque) ;
                    % If #j is periodic
                    if BC.BC{j}.type == 1
                        dfree = BC.BC{j}.ddlfree ;
                        % Replace "free" dof blocked by #j with matching
                        % (actually) free dof of #j
                        d1{k}(id1(match)) = dfree(idfree(match)) ;
                    % If #j is heterogeneous Dirichlet
                    elseif BC.BC{j}.type==0 && BC.BC{j}.ishomogeneous==0
                        % Find corresponding value imposed by #j
                        dval = BC.BC{j}.value ;
                        dval = dval(idfree(match)) ;
                        % Assign it to the matching blocked dof in v
                        d3 = d2 ;
                        d3{k} = BC.BC{j}.ddlbloque ;
                        d3{k} = d3{k}(idfree(match)) ;
                        v(d3{:}) = dval ;
                        % Forget about this blocked dof (now processed)
                        % and the associated "free" dof (blocked by #j)
                        d1{k}(idfree(match)) = [] ;
                        d2{k}(idfree(match)) = [] ;
                    end
                end
            end
            % % END DEBUG % %
            v(d2{:})=v(d1{:});
        elseif BC.BC{i}.type==0 && BC.BC{i}.ishomogeneous==0
            d1{k}=BC.BC{i}.ddlbloque;
            v(d1{:})=BC.BC{i}.value;
        end
    end
    
elseif isa(u,'PCMATRIX')
    if k~=1
        error('')
    end
    s=size(u);
    s(k)=BC.nbddl;
    v=double(u);
    v=unfreevector(BC,v,k);
    v=PCMATRIX(v,s,getPC(u));
elseif isa(u,'PCRADIAL')
    v=funV(u,@(V)unfreevector(BC,V,k));
elseif isa(u,'PCRADIALMATRIX')
    v=setV(u,unfreevector(BC,getV(u),k));
elseif isa(u,'MULTIMATRIX')
    if k~=1
        error('')
    end
    if iscell(u)
        v = setvalue(u,unfreevector(BC,getvalue(u)));
        v = setsize(v,[BC.nbddl,size(u,2)]);
    else
        s=size(u);
        v=reshape(double(u),[s(1),s(2)*length(u)]);
        v=unfreevector(BC,v,k);
        v=reshape(v,[BC.nbddl*s(2),length(u)]);
        u=setvalue(u,v);
        u=setsize(u,[BC.nbddl,size(u,2)]);
        v=u;
    end
elseif isa(u,'FENODEFIELD')
    
    warning('unfreevector peut poser un probleme avec un FENODEFIELD')
    v = setvalue(u,unfreevector(BC,getvalue(u),k));
    
elseif isa(u,'cell')
    v=u;
    for ll=1:length(u)
        v{ll}=unfreevector(BC,u{ll},k) ;
    end
else
    
    v=unfreevector(u,BC,k) ;
    
end


