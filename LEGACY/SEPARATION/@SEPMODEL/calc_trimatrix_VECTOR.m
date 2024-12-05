function A=calc_trimatrix_VECTOR(SM,a,varargin)
% 2 - Sortie de type vecteur :
%     \int{ grad(u) v w } || \int{ u grad(v) w } || \int{ u v grad(w) }
%     Ici, grad(u) sera vectoriel : A={SEP}


p=getp(a);
map=SM.mapping;
spdim =spatialdim(SM);
nspdim=1:SM.dim;
nspdim(spdim)=[];
A=cell(length(spdim),1);

	%-------------------------%
	% Initialisation de Kappa %
    %-------------------------%
    
if isempty(map)
    % Si il n'y a pas de mapping
    K=sepone(ones(1,SM.dim)); % SYMBOLIQUE
else
    if spatial_mapping(map)
        K=getDetFFm1(map);
    else
        K=getDetF(map);
    end
end

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



if isa(K,'SEP')
    % int( k grad(u)vw ) 
    % l'utilisateur doit realiser {A}*{v} correctement.
    OpSp=cell(1,length(spdim));
    for d=1:length(spdim), OpSp{1,d} = getfun(SM,[0 0 0],spdim(d),varargin{:});end
    OpSp=repmat({OpSp(1,:)},length(spdim),1);
    for d=1:length(spdim), OpSp{d}{1,d} = getfun(SM,p,spdim(d),varargin{:});end
    % Finalisation
    for i=1:length(spdim)
        A{i}=finalisation(SEP(OpSp{i}),OpSt,K);
    end
    A=expand(A);
elseif isa(K,'cell') && all(size(K)==[length(spdim) length(spdim)])
    % int( K.grad(u)vw ) 
    % l'utilisateur doit realiser {A}*{v} correctement.
    OpSp=cell(1,length(spdim));
    for d=1:length(spdim), OpSp{1,d} = getfun(SM,[0 0 0],spdim(d),varargin{:});end
    OpSp=repmat({OpSp(1,:)},length(spdim),1);
    for d=1:length(spdim), OpSp{d}{d} = getfun(SM,p,spdim(d),varargin{:});end
    % Finalisation
    for i=1:length(spdim)
        
        % Reperer les K{i,j} non vides
        J=[];
        for j=1:length(spdim)
            if ~(isa(K{i,j},'double') && K{i,j}==0)
                J=[J j];
            end
        end
        % Sommer sur les J
        j=1;
        A{i}=finalisation(SEP(OpSp{J(j)}),OpSt,K{i,J(j)});
        for j=2:length(J)
            A{i}=A{i}+expand(finalisation(SEP(OpSp{J(j)}),OpSt,K{i,J(j)}));
        end
    end
else
    error('Cas non repertorie')
end

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


function Aex=expand(A)
% Dans le cas ou une dimension spatiale est >1,
% il faut reorganiser la structure de A

for i=1:numel(A)
    C=A{i}.F{i};
    if ~isa(C,'cell'), C={C};end
    Ai=cell(size(C));
    for k=1:numel(Ai)
        Ai{k}=A{i};
        Ai{k}.F{i}=C{k};
    end
    A{i}=Ai;
end

Aex={};
for i=1:numel(A)
    Aex=[Aex;A{i}];
end


end




