function x = sumorthogleaves(Q,W)
% function x = sumorthogleaves(Q,W)

if max(rank(W))~=1
    error('max(rank(W)) must be 1')
end


d2i = Q.dim2ind;

x = Q+W;
% Computation of alpha_i^\mu = <W^\mu,Q^\mu>_\mu
alpha = cellfun(@(x,y) full(x'*y),W.U(d2i),Q.U(d2i),'uniformoutput',false);

% Orthogonal component:
% RW^mu = W^mu - \sum_{i=1}^{r_mu} \alpha_i^mu Q_i^\mu
RW = W.U(d2i);
RW = cellfun(@(x,y,z) x-z*y',RW,alpha,Q.U(d2i),'uniformoutput',false);

normRW = cellfun(@(x) norm(x,'fro'),RW,'uniformoutput',false);
RW = cellfun(@(x,y) x/y,RW,normRW,'uniformoutput',false);


in_space = cellfun(@(x,y) x/norm(y,'fro')<4*eps,normRW,W.U(d2i));

q = cell(1,ndims(Q));
q(in_space) = Q.U(d2i(in_space));
q(~in_space) = cellfun(@(x,y) [x y], Q.U(d2i(~in_space)),RW(~in_space),...
    'uniformoutput',false);
r = cellfun(@(x,y) full(x'*y),x.U(d2i),q,'uniformoutput',false);

% x.U = U;
for mu = 1:ndims(Q)
    x.U{d2i(mu)} = q{mu};
end

x_parent = x.parent;
x_is_left = x.is_left;
for mu = 1:ndims(Q)
    ii = d2i(mu);
    ii_par = x_parent(ii);
    if(x_is_left(ii))
        % left child of parent node
        x.B{ii_par} = ttm(x.B{ii_par}, r{mu}', 1);
    else
        % right child of parent node
        x.B{ii_par} = ttm(x.B{ii_par}, r{mu}', 2);
    end
end

end

