function be = eval_elem(a,elem,node,varargin)
% function be = eval_elem(a,elem,node,varargin)

if ischarin('parent',varargin) % compute "be" on parent model instead
   parentModel = getcharin('parent',varargin) ;
   elemChild = elem ;
   elem = getgroupelem(parentModel,getparam(elem,'parentgroupelem')) ;
   node = parentModel.node ;
else % Just in case, yet unlikely
    elemChild = elem ;
    warning('No parent model provided. Will attempt to assemble anyway.')
end

xnode = node(elem);
if ischarin('quadorder',varargin)
    gaussChild = calc_gauss(a,elemChild,getcharin('quadorder',varargin));
else
    gaussChild = calc_gauss(a,elemChild);
end
gauss = change_gauss_syscoord(gaussChild,elemChild,elem) ;
xgauss = gauss.coord;
p = getp(a);
pk = getpk(a);

if any([p,pk]==0)
    N = calc_N(elem,xnode,xgauss);
end
if any([p,pk]==1)
    [DN,~] = calc_B(elem,xnode,xgauss);
    detJ = calc_detJ(elemChild,node(elemChild),gaussChild.coord);
else
    detJ = calc_detJ(elemChild,node(elemChild),gaussChild.coord);
end

if ischarin('V1',varargin)
    V1 = getcharin('V1',varargin) ;
else
    V1 = [] ;
end
if ischarin('V2',varargin)
    V2 = getcharin('V2',varargin) ;
else
    V2 = [] ;
end
V = {V1,V2};

repempty = find([isempty(V1),isempty(V2)]);
% FORMATER [V]
for i=1:2
    if isempty(V{i})
        % V{i}=':'
        if p(i)==0
            V{i} = N;
        elseif p(i)==1
            V{i} = DN;
        end
    else
        % V{i}=v
        if p(i)==0
            if isa(V{i},'cell') && isa(V{i}{1},'function_handle')
                V{i} = V{i}{1}(xglob);
            else    V{i} = N*localize(elem,V{i});
            end
        elseif p(i)==1
            if isa(V{i},'cell') && isa(V{i}{1},'function_handle')
                V{i} = V{i}{1}(xglob);
                V{i} = MYDOUBLEND(V{i});
            else    V{i} = DN*localize(elem,V{i});
            end
        end
    end
end

% FORMATER [K]
k = getk(a);
if ~isempty(k) && ~isempty(pk)
    switch pk
        case 0
            if isa(k,'function_handle')
                cas = 'ELEM';
                if cas=='NODE'
                    error('Non programme')
                elseif cas=='ELEM'
                    error('Non programme')
                end
            elseif isa(k,'cell') && getdim(elem)==length(k)
                % Vecteur k
                dim = getdim(elem);
                nm=size(k);
                for d=1:numel(k)
                    k{d} = mylocnode(k{d},elem,node,N);
                end
                K = k{1};
                for d=2:numel(k)
                    K = [K;k{d}];
                end
                k = K;
                k = reshape(k,[nm,sizeND(K)]);
            else
                k = mylocnode(k,elem,node,N);
            end
        case 1
            k = DN*localize(elem,k);
    end
end

rep1 = find(p==1);
if ~isempty(k) && mod(length(rep1),2)==1
    V{rep1(1)} = k' * V{rep1(1)};
    p(rep1(1)) = 0;
elseif ~isempty(k) && mod(length(rep1),2)==0
    V{1} = k * V{1};
end


% if size(V{1})==size(V{2})
% Si il y en a une, la sommation grad(~).grad(~)
% est realisable. be sera un MYDOUBLEND (format
% classique)
be = V{1}'*V{2};

% fit to child model
if ischarin('parent',varargin)
    
    % List parent elements to keep (because they have children here)
    c2pElemList = find_edges(getconnec(elemChild),getconnec(elem)) ;
    %     c2pElemList = getparam(elemChild,'parentelemnum') ; % wrong
    c2pElemList = [getnumelem(elemChild) c2pElemList] ;
    % Locate parent elements listed previously
    isBorderElem = ismember(getnumelem(elem),c2pElemList(:,2)) ;
    isBorderElem = reshape(MYDOUBLEND(isBorderElem),1,1,[]) ;
    % Set to zero for non-border elements
    be = be*isBorderElem ;
    
    % Store duplicates (nth children with n>1) separately
    c2pElem = [] ;
    dup = (1:size(c2pElemList,1))' ; % potential duplicates
    while ~isempty(dup)
        [~,o2u] = unique(c2pElemList(dup,2),'stable') ; % get original to unique index
        c2pElem = [c2pElem ; {c2pElemList(dup(o2u),:)}] ;
        dup(o2u) = [] ;
    end
    % Read: element with numero c2pElem{n}(k,2) has element with
    % numero c2pElem{n}(k,1) as its nth child (for any k).
    
    % Assign children's jacobian to the right parents
    % (Cf. change_gauss_syscoord)
    detJChild = detJ ;
    szJ = size(detJ) ;
    if size(szJ,2)<4 % if, e.g., only one gauss point per element
        szJ = [szJ ones(1,4-size(szJ,2))];
    end
    nGauss = szJ(4) ;
    szJ(3:4) = [size(be,3) szJ(4)*numel(c2pElem)] ;
    detJ = MYDOUBLEND(zeros(szJ)) ;
    for k = 1:numel(c2pElem)
        iChild = ismember(c2pElemList(:,1),c2pElem{k}(:,1)) ;% same order because 'stable' above
        % Adapt to number of Gauss points per element
        currentPoints = (nGauss*k-(nGauss-1)):(nGauss*k) ;
        detJ(:,:,c2pElem{k}(:,2),currentPoints) = detJChild(:,:,iChild,:) ;
    end
end

be = sum(gauss.w*abs(detJ)*be,4);


if isempty(repempty)
    be = double(sum(be,3));
end

% elseif ischarin('derivate_spdim',varargin)
%     dim2keep = calc_gauss(a,elem,getcharin('quadorder',varargin));
%     % Le cas a concidere ici est le suivant :
%     %   \int grad(u) v :
%     %   b=BILINFORM(1,0); b{M}{:,:}
%     % Dans ce cas, il n'y a qu'un seul gradient.
%     % Si la dimension spatiale est >1, la sommation
%     % n'est pas possible. Le resultat sera alors une
%     % cellule de MYDOUBLEND :
%     %   be={ \int u,x1.v ; \int u,x2.v ... }
%     %           be1            be2     ...
%     sv1 = size(V{1});
%     sv2 = size(V{2});
%     sv = [sv1(1) sv2(1)];
%     if min(sv)~=1
%         error('format non prevu...')
%     end
%
%     [sv,rep] = max(sv);
%     % attention a la permutation V{1} V{2}... rep!
%     if rep~=2
%         Vtmp = V{2};
%         V{2} = V{1};
%         V{1} = Vtmp;
%     end
%     be = cell(sv,1);
%
%     for d=1:sv
%         be{d} = V{1}'*V{2}(d,:,:,:);
%         be{d} = sum(gauss.w*abs(detJ)*be{d},4);
%         if isempty(repempty)
%             be{d} = double(sum(be{d},3));
%         end
%     end
% end





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
    tau = MYDOUBLEND(tau);
else
    tau = tau(getnumber(elem),:);
    tau = reshape(full(tau)',size(tau,2),1,getnbelem(elem));
    tau = MYDOUBLEND(tau);
end

return

