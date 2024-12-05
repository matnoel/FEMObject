function be = eval_elem(a,elem,node,varargin)
% function be = eval_elem(a,elem,node,V1,...,Vn,varargin)

V = varargin(1:getn(a));
varargin = varargin(getn(a)+1:end);

xnode = node(elem);
if ischarin('quadorder',varargin)
    gauss = calc_gauss(a,elem,getcharin('quadorder',varargin));
else
    gauss = calc_gauss(a,elem);
end
xgauss = gauss.coord;
xglob = calc_x(elem,xnode,xgauss);

if any([a.p,a.pk]==0)
    N = calc_N(elem,xnode,xgauss);
end
if any([a.p,a.pk]==1)
    [DN,detJ] = calc_DN(elem,xnode,xgauss);
else
    detJ = calc_detJ(elem,xnode,xgauss);
end

% FORMATER [:]
repempty = find(isemptyarg(V));
repnonempty = find(~isemptyarg(V));
for i=1:length(repempty)
    j = repempty(i);
    switch a.p(j)
        case 0
            V{j} = N;
        case 1
            V{j} = DN;
    end
end

% FORMATER [V]
for i=1:length(repnonempty)
    j = repnonempty(i);
    
    switch a.p(j)
        case 0
            if isa(V{j},'cell') && isa(V{j}{1},'function_handle')
                V{j} = V{j}{1}(xglob);
            else
                V{j} = N*localize(elem,V{j});
            end
        case 1
            if isa(V{j},'cell') && isa(V{j}{1},'function_handle')
                V{j} = V{j}{1}(xglob);
                V{j} = MYDOUBLEND(V{j});
            else
                V{j} = DN*localize(elem,V{j});
            end
        otherwise
            error('pas prevu')
    end
end

% FORMATER KAPPA : Soit champ par noeud, elem ou bien function_handle
if ~isempty(a.k) && ~isempty(a.pk)
    switch a.pk
        case 0
            if isa(a.k,'function_handle')
                cas = 'ELEM';% ou 'NODE'
                if cas=='NODE'
                    C = getcoord(node);
                    EVAl = zeros(size(C,1),1);
                    for i=1:size(C,1)
                        EVAl(i) = a.k(C(i,:));
                    end
                    a.k = mylocnode(EVAl,elem,node,N);
                elseif cas=='ELEM'
                    C = getcoord(node);
                    EVAL = zeros(getnbelem(elem),1);
                    for i=getnumber(elem)' % ATTENTION : peut etre erreurs
                        eval_coord_i = mean(C(getnumnode(getelem(elem,i)),:),1);
                        EVAL(i) = a.k(eval_coord_i);
                    end
                    a.k = mylocelem(EVAL,elem);
                end
            elseif isa(a.k,'cell') && getdim(elem)==length(a.k)
                % Vecteur k
                
                nm=size(a.k);
                for d=1:numel(a.k)
                    try
                    a.k{d} = mylocnode(a.k{d},elem,node,N);
                    catch
                    a.k{d} = mylocelem(a.k{d},elem);
                    end
                end
                K = a.k{1};
                for d=2:numel(a.k)
                    K = [K;a.k{d}];
                end
                a.k = K;
                a.k = reshape(a.k,[nm,sizeND(K)]);
            else
                try
                    a.k = mylocnode(a.k,elem,node,N);
                catch
                    a.k = mylocelem(a.k,elem);
                end
            end
            
        case 1
            if isa(a.k,'function_handle')
                % RESTE A FAIRE...
            else
                a.k = DN*localize(elem,a.k);
            end
    end
end

rep1 = find(a.p==1);
if ~isempty(a.k) && mod(length(rep1),2)==1
    V{rep1(1)} = a.k' * V{rep1(1)};
    if length(repempty)~=3
        a.p(rep1(1)) = 0;
    end
elseif ~isempty(a.k) && mod(length(rep1),2)==0
    V{1} = a.k * V{1};
end


if length(repempty)==0
    % cas scalaire
    be = prodV(V,a.p);
elseif length(repempty)==1
    % cas vecteur
    repl = 1:repempty(1)-1;
    repr = repempty(1)+1:length(V);
    [bel,gradl] = prodV(V(repl),a.p(repl));
    [ber,gradr] = prodV(V(repr),a.p(repr));
    if a.p(repempty(1))==0
        be = V{repempty(1)}'*(bel'*ber);
    else
        if gradr && ~gradl
            be = bel*V{repempty(1)}'*ber;
        elseif gradl && ~gradr
            be = (bel'*V{repempty(1)})'*ber;
        elseif ~gradl && ~gradr
            be = bel*V{repempty(1)}'*ber;
        end
    end
elseif  length(repempty)==2
    % cas matrice
    repl = 1:repempty(1)-1;
    repm = repempty(1)+1:repempty(2)-1;
    repr = repempty(2)+1:length(V);
    [bel,gradl] = prodV(V(repl),a.p(repl));
    [bem,gradm] = prodV(V(repm),a.p(repm));
    try
        [ber,gradr] = prodV(V(repr),a.p(repr));
    catch
        keyboard
    end
    if a.p(repempty(1))==0 && a.p(repempty(2))==0
        if ~gradl
            be = (V{repempty(1)}'*V{repempty(2)})*(bel*(bem'*ber));
        else
            be = (V{repempty(1)}'*V{repempty(2)})*(bel'*(bem*ber));
        end
    elseif a.p(repempty(1))==1 && a.p(repempty(2))==0
        
        if gradl
            be1 = (bel'*V{repempty(1)})'*V{repempty(2)};
            be2 = bem'*ber;
        elseif gradm
            be1 = bel*(V{repempty(1)}'*bem);
            be2 = V{repempty(2)}*ber;
        elseif gradr
            be1 = V{repempty(1)}'*ber;
            be2 = bel*V{repempty(2)}*bem;
        else
            be1 = bel*V{repempty(1)}'*bem;
            be2 = V{repempty(2)}*ber;
        end
        be= be1*be2;
    elseif a.p(repempty(1))==0 && a.p(repempty(2))==1
        if gradr
            be2 = (ber'*V{repempty(2)});
            be1 = (bel'*bem)*V{repempty(1)}';
        elseif gradm
            be2 = (bem'*V{repempty(2)})*bel;
            be1 = bel*V{repempty(1)}';
        elseif gradl
            be1 = V{repempty(1)}';
            be2 = (bel'*V{repempty(2)})*bem*ber;
        else
            be1 = bel*V{repempty(1)}'*bem;
            be2 = V{repempty(2)}*ber;
        end
        be = be1*be2;
    elseif a.p(repempty(1))==1 && a.p(repempty(2))==1
        if gradl
            be1 = (bel'*V{repempty(1)})';
            if gradm
                be2 = (bem'*V{repempty(2)})*ber;
            elseif gradr
                be2 = bem*(ber'*V{repempty(2)});
            end
            be = be1*be2;
        elseif gradm
            be1 = bel*(V{repempty(1)}'*bem);
            be2 = ber'*V{repempty(1)};
            be = be1*be2;
        else
            be = bel*(V{repempty(1)}'*bem*V{repempty(2)})*ber;
        end
    end
elseif  length(repempty)==3
    % cas tenseur d'ordre 3
    % Retourne le MYDOUBLEND be{ n1 n2 n3 elem gauss}
    % But : stocker cette info dans un SPARSETENSOR afin d'eviter la
    % repetition de trop de calculs lourds.
    
    % AU PLUS SIMPLE : cas TRILINFORM UNIQUEMENT
    % cas matrice
    p = getp(a);
    pk = getpk(a); if isempty(pk), pk = 0;end
    rep1 = find(p==1);
    rep0 = find(p==0);
    dim = getdim(elem);
    
    if (isempty(rep1) && pk==0) || (sum([p pk])==1 && dim==1)
        % kuvw OU kuv(w,x)
        be = gauss.w*abs(detJ)*V{1}'*V{2};
        be = reshape(be,[ 1 size(V{1},2)^2 sizeND(V{1}) ]);
        be = be'*V{3};
        be = reshape(be,[size(V{3},2) size(V{1},2) size(V{1},2) sizeND(V{1})]);
    elseif length(rep1)==1 && pk==1 && ~ischarin('variable2derivate',varargin)
        % grad(k)grad(u)vw , a une permutation (uvw) pres
        be = gauss.w*abs(detJ)*(V{rep1(1)})'*V{rep0(1)};
        be = reshape(be,[ 1 size(V{1},2)^2 sizeND(V{1}) ]);
        be = be'*V{rep0(2)};
        be = reshape(be,[size(V{rep0(2)},2) size(V{rep0(2)},2) size(V{rep0(2)},2) sizeND(V{rep0(2)})]);
        % Attention a la permutation:
        be = ipermute(be,[rep1(1) rep0(1) rep0(2) 4:5]);
    elseif length(rep1)==2 && pk==0 && ~ischarin('variable2derivate',varargin)
        % kgrad(u)grad(v)w , a une permutation (uvw) pres
        be = gauss.w*abs(detJ)*V{rep1(1)}'*V{rep1(2)};
        be = reshape(be ,[ 1 size(V{1},2)^2 sizeND(V{1}) ]);
        be = be'*V{rep0(1)};
        be = reshape(be,[size(V{1},2) size(V{1},2) size(V{1},2) sizeND(V{1})]);
        % Attention a la permutation :
        be = ipermute(be,[rep1 rep0 4:5]);
    elseif length(rep1)==1 && dim~=1 && ischarin('variable2derivate',varargin)
        % kuv(w,xi) ET dim>2  , a une permutation (kuvw) pres
        % La variable de derivation est :
        var2der = getcharin('variable2derivate',varargin);
        V{rep1} = V{rep1}(var2der,:,:,:);
        
        % finalisation, exactement comme le cas (kuvw OU kuv(w,x))
        be = gauss.w*abs(detJ)*V{1}'*V{2};
        be = reshape(be,[ 1 size(V{1},2)^2 sizeND(V{1}) ]);
        be = be'*V{3};
        be = reshape(be,[size(V{3},2) size(V{1},2) size(V{1},2) sizeND(V{1})]);
    elseif length(rep1)==2 && dim~=1 && ischarin('variable2derivate',varargin)
        % ku(v,xi)(w,xj) ET dim>2  , a une permutation (kuvw) pres
        % La variable de derivation est :
        var2der = getcharin('variable2derivate',varargin);
        V{rep1(1)} = V{rep1(1)}(var2der(1),:,:,:);
        V{rep1(2)} = V{rep1(2)}(var2der(2),:,:,:);
        
        % finalisation, exactement comme le cas (kuvw OU kuv(w,x))
        be = gauss.w*abs(detJ)*V{1}'*V{2};
        be = reshape(be,[ 1 size(V{1},2)^2 sizeND(V{1}) ]);
        be = be'*V{3};
        be = reshape(be,[size(V{3},2) size(V{1},2) size(V{1},2) sizeND(V{1})]);
    else
        error('Cas non repertorie')
    end
    
    be = sum(be,5);
    return
end


be = sum(gauss.w*abs(detJ)*be,4);

% Cas scalaire : finir le travail de sommation
if isempty(repempty)
    be = double(sum(be,3));
end



function [be,grad] = prodV(V,p)
% function [be,grad] = prodV(V,p)

if isempty(V)
    be = 1;
    grad = false;
else
    be = V{1};
    grad = (p(1)==1);
    for i=2:length(V)
        if p(i)==1 && grad
            be = be'*V{i};
            grad = false;
        elseif p(i)==1 && ~grad
            be = be*V{i};
            grad = true;
        elseif p(i)==0
            be = be*V{i};
        end
    end
end

return

% 
% function c = mylocnode(c,elem,N)
% % function c = mylocnode(c,elem,N)
% 
% c = full(c);
% if isempty(c)
%     c = 1;
% elseif numel(c)~=getdim(elem)
%     c = N*localize(elem,c);
% end
% return
% 



function c = mylocnode(c,elem,node,N)
% function c = mylocnode(c,elem,node,N)

c = full(c);
if isempty(c)
    c = 1;
elseif numel(c)~=getdim(elem)
    if size(c,1)==getnbnode(node)
        c = N*localize(elem,c);
    else
        c = mylocelem(c,elem);
        if size(c,1)==getdim(elem)^2
            c = reshape(c,[getdim(elem),getdim(elem),sizeND(c)]);
        end
    end
    if size(c,1)==1
        c = c';
    end
end
return


function tau = mylocelem(tau,elem)
% function tau = mylocelem(tau,elem)

if numel(tau)==1
    tau = MYDOUBLEND(full(tau));
else
    tau = tau(getnumber(elem));
    tau = reshape(full(tau),1,1,getnbelem(elem));
    tau = MYDOUBLEND(tau);
end

return


