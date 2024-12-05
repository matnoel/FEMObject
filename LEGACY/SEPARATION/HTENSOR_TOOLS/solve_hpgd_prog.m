function [u,result] = solve_hpgd_prog(A,b,solver,varargin)
% function [u,result] = solve_hpgd_prog(A,b,solver,varargin)
%   Progressive HPGD

dim = ndims(A);

param = getparam(solver);

if isfield(param,'timer')
    if param.timer ;
        timer = 1;
        result.time = zeros(1,param.maxorder);
        tic
    else
        timer = 0;
    end
else
    timer = 0;
end

if ~isfield(param,'nodemaxiter')
    param.nodemaxiter = 10;
end
if ~isfield(param,'nodetol')
    param.nodetol = 1e-6;
end


sz = size(b);
b=full_leaves(b);
bu = b;

u=htensor(sz,A.children,A.dim2ind);

erriter = zeros(1,param.maxiter);
errorder = zeros(1,param.maxorder);

normref = norm(b);

dim2ind = A.dim2ind;

for i = 1:param.maxorder
    fprintf('order #%d \n',i)
    % U0 = htenrandn(sz,'orthog',ones(1,2*dim-1),A.children,A.dim2ind);
    U0 = htenones(sz,A.children,A.dim2ind);
    U = U0;
    alpha0 = norm(U0);
    AU = apply_mat_to_vec(A,U);
    UAU = apply_mat_to_mat(U,AU,sz);
    Ubu = apply_mat_to_vec(U,bu);
    
    for kkk = 1:param.maxiter
        for j = 1:dim
            fprintf('%d ',j);
            jj = dim2ind(j);
            
            UAU.U{jj} = A.U{jj};
            % M = reshape(full(UAU),sz(j),sz(j));
            % Previous line is OK, but it is memory consuming
            M = reduce(UAU,j);
            M = reshape(M,sz(j),sz(j));
            
            Ubu.U{jj} = bu.U{jj};
            % S = reshape(full(Ubu),sz(j),1);
            S = reduce(Ubu,j);
            
            U.U{jj} = M\S;
            alpha = norm(U.U{jj});
            U.U{jj} = (1/alpha) * U.U{jj};
            
            AU.U{jj} = prod_one_leaf_mv(A.U{jj},U.U{jj});
            UAU.U{jj} = prod_one_leaf_mm(U.U{jj},AU.U{jj},sz(j));
            Ubu.U{jj} = prod_one_leaf_mm(U.U{jj},bu.U{jj},sz(j));
            
        end
        fprintf('\n');
        erriter(kkk)=abs(alpha-alpha0)/alpha0;
        result.alpha{i}(kkk)=alpha;
        
        fprintf('  iteration #%d - stagnation = %1.3d\n',kkk,erriter(kkk))
        if erriter(kkk) < param.itercrit
            break
        end
        
        alpha0 = alpha;
    end
    
    if i==1
        u = orthog(alpha*U);
    else
        u = orthog(u+alpha*U);
        optsu.max_rank = max(rank(u));
        optsu.rel_eps = 1e-10;
        
        fprintf('  NODES Update\n')
        for k=1:param.nodemaxiter
            u0 = u;
            for t=u.nr_nodes:-1:1
                if ~u.is_leaf(t) || param.updateleaves
                    % if t~=4
                    fprintf('    Node #%d\n',t)
                    % u = truncate(u,optsu);
                    u = orthog(u);
                    u = solve_updatenode(A,b,u,t);
                    % end
                end
            end
            fit = norm(orthog(u-u0))/norm(u0);
            fprintf('     Stagnation = %d\n',fit);
            if fit < param.nodetol
                break
            end
        end
    end
    optsu.max_rank = max(rank(u));
    optsu.rel_eps = 1e-10;
    %u = truncate(u,optsu);
    u = orthog(u);
    bu = orthog(b-apply_mat_to_vec(A,u));
    errorder(i) = norm(bu)/normref;
    fprintf(' order #%d - error = %d\n',i,errorder(i));
    if param.storeiter==1
        result.u{i}=u;
    end
    if timer
        result.time(i) = toc;
    end
end

result.error = errorder;

end


function Ax = prod_one_leaf_mv(A,x)
k_A = size(A, 2);
sz = size(x,1);
AxU = cell(1,k_A);
for jj=1:k_A
    Ajj = reshape(A(:, jj), sz, sz);
    AxU{jj} = Ajj * x;
end
Ax = cell2mat(AxU);
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

