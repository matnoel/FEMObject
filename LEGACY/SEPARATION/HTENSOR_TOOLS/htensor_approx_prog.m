function [u,result] = htensor_approx_prog(b,varargin)
% function [u,result] = htensor_approx_prog(b,varargin)

dim = ndims(b);
d2i = b.dim2ind;


solver = SEPSOLVER(dim,varargin{:});
param = getparam(solver);

n = size(b);
bu = b;
nnodes = 2*dim-1;

erriter = zeros(1,param.maxiter);
errorder = zeros(1,param.maxorder);

normref = norm(b);

I = sepopereye(n);
I = sepmatrixtohtensor(I);

optsu.rel_eps = 1e-12;

for i=1:param.maxorder
    if param.display
        fprintf('order #%d \n',i)
    end
    
    U0 = htenrandn(n,'orthog',ones(1,nnodes));
    alpha0 = 1;
    U = U0;
    Ubu = apply_mat_to_vec(U,bu);
    
    for kkk=1:param.maxiter
        for j=1:dim
            jj = d2i(j);
            Ubu.U{jj} = bu.U{jj};
            S = full(Ubu);
            alpha = norm(S(:));
            U.U{jj} = S(:)./alpha;
            Ubu.U{jj} = prod_one_leaf_mm(U.U{jj},bu.U{jj},n(j));
        end
        
        erriter(kkk)=abs(alpha-alpha0)/(alpha+alpha0);
        if param.display
            fprintf('  iteration #%d - stagnation = %d\n',kkk,erriter(kkk))
        end
        
        if kkk > 1 && (erriter(kkk)<param.itercrit)
            break
        end
        alpha0=alpha;
    end
    
    if i==1
        u = U;
    else
        u = orthog(u + alpha*U);
    end
    
    optsu.max_rank = max(rank(u));
    for kkk = 1:param.nodemaxiter
        fprintf('  NODES Update\n')
        u0 = u;
        for t=u.nr_nodes:-1:1
            if ~u.is_leaf(t) || param.updateleaves
                fprintf('    Node #%d\n',t)
                u = truncate(u,optsu);
                u = orthog(u);
                u = solve_updatenode(I,b,u,t);
            end
        end
        fit = norm(orthog(u-u0))/norm(u0);
        fprintf('     Stagnation = %d\n', fit);
        if fit < param.nodetol
            break
        end
    end
    bu = orthog(b-u);
    errorder(i)=norm(bu)/normref;
    
    if param.display
        fprintf('  order #%d - error = %d \n',i,errorder(i))
    end
    if errorder(i) < param.tol
        break
    end
end

result.errorder = errorder;

end

function Ax = prod_one_leaf_mm(A,x,sz)
k_A = size(A, 2);
AxU = cell(1, k_A);
for jj=1:k_A
    Ajj = reshape(A(:, jj), 1, sz);
    AxU{jj} = Ajj * x;
end
Ax = cell2mat(AxU);
end
