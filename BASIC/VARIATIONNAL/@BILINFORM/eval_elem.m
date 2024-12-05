function be = eval_elem(a,elem,node,V1,V2,varargin)
% function be = eval_elem(a,elem,node,V1,V2,varargin)

xnode = node(elem);
if ischarin('quadorder',varargin)
    gauss = calc_gauss(a,elem,getcharin('quadorder',varargin));
else
    gauss = calc_gauss(a,elem);
end
xgauss = gauss.coord;
p = getp(a);
pk = getpk(a);

if any([p,pk]==0)
    N = calc_N(elem,xnode,xgauss);
end
if any([p,pk]==1)
    [DN,detJ] = calc_B(elem,xnode,xgauss);
else
    detJ = calc_detJ(elem,xnode,xgauss);
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
% k = getk(a);
% if ~isempty(k) && ~isempty(pk)
%     switch pk
%         case 0
%             if isa(k,'cell') && isa(k{1},'function_handle')
%                 k = k{1}(xglob);
%             else
%                 k = N*localize(elem,k);
%             end
%         case 1
%             if isa(k,'cell') && isa(k{1},'function_handle')
%                 k = k{2}(xglob);
%                 k = MYDOUBLEND(k);
%             else
%                 k = DN*localize(elem,k);
%             end
%         otherwise
%             error('pas prevu')
%     end
% end
% k = getk(a);
% if ~isempty(k) && ~isempty(pk)
%     switch pk
%         case 0
%             if isa(k,'function_handle')
%                 K = zeros(size(V{1},1),size(V{1},1),size(V{1},3),size(V{1},4));
%                 XG = N*xnode;
%                 for i=1:size(V{1},3)  % Sur les elements
%                     XGe      = zeros(size(XG,2),size(XG,4));
%                     XGe(:,:) = XG(1,:,i,:) ;
%                     K(:,:,i,:) = k(XGe');
%                 end
%                 k = MYDOUBLEND(K);
%             else
%                 k = mylocnode(k,elem,node,N);
%             end
%         case 1
%             k = DN*localize(elem,k);
%     end
% end



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

