function A=calc_trimatrix_SCALAR(SM,a,varargin)
% 1 - Sortie de type scalaire :
%     \int{ u v w }
%     \int{ grad(u).grad(v) w }  || \int{ u grad(v).grad(w) }
%     Les cas 2 avec une dimension spatiale=1 sont aussi a prendre en
%     compte.
%     La sortie est alors une (H)SEP


p=getp(a);
map=SM.mapping;
spdim =spatialdim(SM);
nspdim=1:SM.dim;
nspdim(spdim)=[];



	%-------------------------%
	% Initialisation de Kappa %
    %-------------------------%
if isempty(map)
    % Si il n'y a pas de mapping
    K=sepone(ones(1,SM.dim)); % SYMBOLIQUE
else
    if spatial_mapping(map)
        if sum(p)==0
            K=getDetF(map);
        elseif sum(p)==1 % Forcement en dimension spatiale=1
            K=getDetFFm1(map);
            K=K{1};
        elseif sum(p)==2
            K=getDetFFmtFm1(map);
        else
            error('Cas non repertorie')
        end
    else
        % Cas stochastique
        K=getDetF(map,SM.dim);
    end
end

% K=getk(a);


    %---------------------------------------------------------%
    % STOCHASTIQUE : plus facile car getmasse a tous le coups %
    %---------------------------------------------------------%
OpSt=cell(1,length(nspdim));
for d=1:length(nspdim)
    %%% Long et lourd si KAPPA est stochastique :
    %OpSt{d}= getfun(SM,[0 0 0],nspdim(d),varargin{:});
    % Rapide si deja calcule :
    OpSt{d}=@(KAPPA) SM.F{nspdim(d)}.masse;
end
OpSt=SEP(OpSt);


    %------------------------------------------------%
	% SPATIAL : Attention aux derivees et au mapping %
    %------------------------------------------------%

if sum(p)==0
    % int( kuww )
    OpSp=cell(1,length(spdim));
    for d=1:length(spdim)
        OpSp{d} = getfun(SM,p,spdim(d),varargin{:});
    end
    OpSp=SEP(OpSp);
    
elseif sum(p)==1 % Forcement en dimension spatiale=1
    % int( kuww,x )
    if length(spdim)~=1, error('Incoherence : Forcement en dimension spatiale=1'); end
    OpSp = getfun(SM,p,[],spdim,varargin{:});
    OpSp = SEP(OpSp);
elseif sum(p)==2
    if isa(K,'SEP')
        % int( k grad(u).grad(v)w )
        OpSp=cell(1,length(spdim));
        for d=1:length(spdim), OpSp{1,d} = getfun(SM,[0 0 0],spdim(d),varargin{:}); end
        OpSp=repmat(OpSp,length(spdim),1);
        for d=1:length(spdim), OpSp{d,d} = getfun(SM,p,spdim(d),varargin{:}); end
        OpSp = SEP(OpSp);
    else
        % int( grad(u).K.grad(v)w )
        P=myindex(p,length(spdim));
        OpSp=cell(length(spdim));
        Sij=[];
        for i=1:length(spdim)
        for j=1:length(spdim)
            % 2-Pour ij fixe, calcul de OpSpij=@(k)(\int u,xi k v,xj w)
            OpSp{i,j}=cell(1,length(spdim));
            for d=1:length(spdim), OpSp{i,j}{d} = getfun(SM,P{i,j,d},spdim(d),varargin{:}); end
            % 3-Finalisation du terme {i,j}
            if (isa(K{i,j},'double') && K{i,j}==0)
                % Si terme nul...
                OpSp{i,j}=[];
            else
                % Si terme existant :
                Sij=[Sij, sub2ind([length(spdim) length(spdim)],i,j)];
                OpSp{i,j}=finalisation(SEP(OpSp{i,j}),OpSt,K{i,j});
            end
        end
        end
        
        %   4-Sommer A = \sum_ij OpSp{i,j}
        
        A=OpSp{Sij(1)};
        if A.dim==2 && all(cellfun(@(o) o.F{1,2}.m,OpSp(Sij(:)))==1)
            for i=Sij(2:end)
                A.F{1,1}= {A.F{1,1} + OpSp{i}.F{1,1}};
            end
        else
            for i=Sij(2:end)
                A=A+OpSp{i};
            end
        end
        
        return
    end
end


    %-------------------------%
    % Finalisation du travail %
    %-------------------------%
A=finalisation(OpSp,OpSt,K);

end


function fun=getfun(SM,p,d,varargin)
fun=@(kappa) calc_trimatrix(MULTILINFORM(p,kappa,0),SM.F{d}.model,varargin{:});
end

function A=finalisation(OpSp,OpSt,K)
if OpSt.dim~=0 && OpSp.dim~=0
    if OpSp.m==1 && OpSt.m==1, A=SEP([OpSp.F OpSt.F]); 
    else                       A=HSEP({OpSp   OpSt  });
    end
    
    % Compatibilite de A et K
    if isa(K,'HSEP')    , A=HSEP(A,K.tree);
    elseif isa(A,'HSEP'), K=change_tree(K,A);%K=HSEP(K,A.tree);
    end
        
    A= mtimeslike(A,K,@(fun,k) fun(k),'');
elseif OpSt.dim==0
    A= mtimeslike(OpSp,K,@(fun,k) fun(k),'');
elseif OpSp.dim==0
    A= mtimeslike(OpSt,K,@(fun,k) fun(k),'');
end
end


function P=myindex(p,dim)
% Pas franchement joli a regarder....

    % 1- IND = [ \delta_id \delta_jd ]
    IND=zeros(dim,dim,dim,2);
    for i=1:dim
    for j=1:dim
    for d=1:dim
        if i==d, IND(i,j,d,1)=1;end
        if j==d, IND(i,j,d,2)=1;end
    end
    end
    end
    
    % 2- [p pk]=[1      0 1      0]    <- exemple
    %  ->[P PK]=[IND(1) 0 IND(2) 0]
    P=repmat({zeros(1,3)},[dim dim dim]);
    locD=find(p==1);
    for i=1:dim,for j=1:dim,for d=1:dim,
    P{i,j,d}(locD(1))=IND(i,j,d,1);
    end, end, end
    for i=1:dim,for j=1:dim,for d=1:dim,
    P{i,j,d}(locD(2))=IND(i,j,d,2);
    end, end, end

end




