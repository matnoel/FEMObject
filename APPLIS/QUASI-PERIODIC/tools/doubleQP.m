function tf = doubleQP(t)
% tf = doubleQP(t)
% Converts TuckerLikeTensor into double following QP format interpretation.

assert(isa(t,'AlgebraicTensor'),'Input must be AlgebraicTensor.')

tcore = double(t.core) ;
tsp = t.space.spaces ;
inds = find(tcore) ;
if isempty(inds) % Safety
    tsz = [prod(t.space.sz') 1] ;
    tf = sparse(tsz(1),tsz(2)) ;
    return
end
subs = ind2NDsub(t.space.dim,inds) ;
isOperator = size(t.sz,1) == 2 ;
if isOperator
    tf = sparse([],[],[],prod(t.sz(1,:)),prod(t.sz(2,:))) ;
else % vector
    tf = sparse([],[],[],prod(t.sz(1,:)),1) ;
end
if isempty(inds) ; return ; end % safety
switch t.order
    case 2
        if isOperator
            for l = 1:size(subs,1)
                tf = tf + tcore(inds(l))*kron(tsp{1}{subs(l,1)},...
                    tsp{2}{subs(l,2)}) ;
            end
        else
            for l = 1:size(subs,1)
                tf = tf + tcore(inds(l))*kron(tsp{1}(:,subs(l,1)),...
                    tsp{2}(:,subs(l,2))) ;
            end
        end
    case 3
        if isOperator
            for y = 1:t.space.dim(1,end)
                M = sparse([],[],[],prod(t.sz(1,1:2)),prod(t.sz(2,1:2))) ;
                for l = find(subs(:,end)==y)'
                    M = M + tcore(inds(l))*kron(tsp{2}{subs(l,2)},tsp{1}{subs(l,1)}) ;
                end
                tf = tf + kron(M,tsp{3}{y}) ;
            end
        else
            for y = 1:t.space.dim(1,end)
                M = sparse([],[],[],prod(t.sz(1,1:2)),1) ;
                for l = find(subs(:,end)==y)'
                    M = M + tcore(inds(l))*kron(tsp{2}(:,subs(l,2)),tsp{1}(:,subs(l,1))) ;
                end
                tf = tf + kron(M,tsp{3}(:,y)) ;
            end
        end
end