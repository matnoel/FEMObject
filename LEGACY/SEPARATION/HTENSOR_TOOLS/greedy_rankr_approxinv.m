function [Q,eQ] = constr_precond(A,opts)
% function [Q,eQ] = constr_precond(A,opts)

Q = cell(opts.maxrank,1);
eQ = zeros(opts.maxrank,1);

n = sqrt(size(A));
I = sepmatrixtohtensor(sepopereye(n));
nI = norm(I);


fprintf('Q : order 1\n');
Q{1} = solve_rank1approxinv(A,opts);
R = full_leaves(I-apply_mat_to_mat(Q{1},A,n));
R = orthog(R);
eQ(1) = norm(R)/nI;
fprintf('Q : order %d -- error %d\n',1,eQ(1));


for i = 2:opts.maxrank
    fprintf('Q : order %d\n',i);
    W = solve_rank1approxinv(A,Q{i-1},opts);
    Q{i} = Q{i-1}+W;
    
    R = I-apply_mat_to_mat(Q{i},A,n);
    eQ(i) = norm(R)/nI;
    fprintf('Q : order %d -- error %d\n',i,eQ(i));
end

end
