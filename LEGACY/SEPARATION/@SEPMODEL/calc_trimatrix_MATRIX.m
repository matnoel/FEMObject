function A = calc_trimatrix_MATRIX(SM,a,varargin)
% function A = calc_trimatrix_MATRIX(SM,a,varargin)
% 2 - Sortie de type matrice :
%     \int{ grad(u) V grad(w) }
%     Avec V une matrice (RdC, eq. du Kappa)
%     L'utilisateur doit realiser {u}*{V}*{w} correctement.

p=getp(a);
map=SM.mapping;
spdim =spatialdim(SM);
nspdim=1:SM.dim;
nspdim(spdim)=[];
A=cell(length(spdim));


%-------------------------%
% Initialisation de Kappa %
%-------------------------%

if isempty(map)
    % Si il n'y a pas de mapping
    K=sepone(ones(1,SM.dim)); % SYMBOLIQUE
else
    if spatial_mapping(map)
        % \int{ grad(u).F^{-t}.V.F^{-1}.grad(w) detF}
        % ICI K est d'ordre 4: K{i,a,b,j} tq :
        %   \int{ u,a V_ij w_b K_aijb} = \int{ grad(u).F^{-t}.V.F^{-1}.grad(w) detF}
        %   Soit K_aijb = F^{-t}_ai F^{-1}_jb det(F)
        K=getDetFFmtFm1_T4(map);
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
    % OpSt{d}= getfun(SM,[0 0 0],nspdim(d),varargin{:});
    % Rapide si deja calcule :
    OpSt{d}=@(KAPPA) SM.F{nspdim(d)}.masse;
end
OpSt=SEP(OpSt);




%------------------------------------------------%
% SPATIAL : Attention aux derivees et au mapping %
%------------------------------------------------%
P=myindex(p,length(spdim));

if isa(K,'SEP')
    % int( k grad(u)Vgrad(w) )
    for i=1:length(spdim)
        for j=1:length(spdim)
            OpSp=cell(1,length(spdim));
            for d=1:length(spdim), OpSp{1,d} = getfun(SM,P{i,j,d},spdim(d),varargin{:});end
            A{i,j}=finalisation(SEP(OpSp),OpSt,K);
        end
    end
    A=expand(A);
else
    % int( grad(u)K.V.Kgrad(w) )
    for i=1:length(spdim)
        for j=1:length(spdim)
            % Remplissage de A{i,j} :
            %
            OpSp=cell(1,length(spdim));
            Aij=cell(length(spdim));
            SAB=[];
            for a=1:length(spdim)
                for b=1:length(spdim)
                    for d=1:length(spdim), OpSp{1,d} = getfun(SM,P{a,b,d},spdim(d),varargin{:});end
                    % Possibilite de presence de termes nuls dans K{a,i,j,b} :
                    if isa(K{a,i,j,b},'double') && K{a,i,j,b}==0
                        Aij{a,b}=[];
                    else
                        % recuperer les indices nons nuls
                        SAB=[SAB, sub2ind([length(spdim) length(spdim)],a,b)];
                        % Finaliser
                        Aij{a,b}=finalisation(SEP(OpSp),OpSt,K{a,i,j,b});
                    end
                end
            end
            Aij=expand(Aij);
            
            A{i,j}=Aij{SAB(1)};
            for sab=SAB(2:end)
                A{i,j}=A{i,j}+Aij{sab};
            end
        end
    end
    
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
    elseif isa(A,'HSEP'), K=change_tree(K,A);% K=HSEP(K,A.tree);
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

for i=1:size(A,1)
    for j=1:size(A,1)
        if i==j
            % diagonale
            C=A{i,i}.F{i};
            if ~isa(C,'cell'), C={C};end
            Aij=cell(size(C));
            for k=1:numel(Aij)
                Aij{k}=A{i,i};
                Aij{k}.F{i}=C{k};
            end
        else
            C1=A{i,j}.F{i};
            C2=A{i,j}.F{j};
            if ~isa(C1,'cell'), C1={C1};end
            if ~isa(C2,'cell'), C2={C2};end
            
            Aij=cell(numel(C1),numel(C2));
            for ie=1:numel(C1)
                for je=1:numel(C2)
                    Aij{ie,je}=A{i,j};
                    Aij{ie,je}.F{i}=C1{ie};
                    Aij{ie,je}.F{j}=C2{je};
                end
            end
        end
        A{i,j}=Aij;
    end
end

Aex={};
for i=1:size(A,1)
    Ai={};
    for j=1:size(A,1)
        Ai=[Ai A{i,j}];
    end
    Aex=[Aex;Ai];
end

end





