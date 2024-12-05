function [Q,eQ] = constr_hprecond(A,opts)
% function [Q,eQ] = constr_hprecond(A,opts)

Q = cell(opts.maxrank,1);
eQ = zeros(opts.maxrank,1);

n = sqrt(size(A));
I = sepmatrixtohtensor(sepopereye(n));
nI = norm(I);

T = 1:A.nr_nodes;
T(A.is_leaf) = [];

warnid = 'MATLAB:nearlySingularMatrix';
warning('off',warnid);



fprintf('Q : order 1\n');

Err = @(x) I-apply_mat_to_mat(x,A,n);

Q{1} = solve_rank1approxinv(A,opts);
Q{1} = orthog(Q{1});
eQ(1) = norm(Err(Q{1}))/nI;
fprintf('Q : order %d -- error %d\n',1,eQ(1));


for i = 2:opts.maxrank
    fprintf('Q : order %d\n',i);
    W = solve_rank1approxinv(A,Q{i-1},opts);
    QQ = sumorthogleaves(Q{i-1},W);
    QQ = orthogtransfertensors(QQ);
    for j = 1:5
        Q0 = QQ;
        for t = T
            fprintf('%d ',t);
            QQ = solve_updatenonleafnodeinv(A,QQ,t,opts);
            QQ = orthogparents(QQ,t);
        end
        fprintf('\n');
        err = norm(QQ-Q0)/norm(QQ);
        fprintf('Node update : %d\n',err);
        if  err < 1e-6
            break
        end
    end
    tr.max_rank = max(rank(QQ));
    tr.rel_eps = 2*eps;
    QQ = truncate(QQ,tr);
    
    Q{i} = QQ;
    
    eQ(i) = norm(Err(Q{i}))/nI;
    fprintf('Q : order %d -- error %d\n',i,eQ(i));
end
warning('on',warnid);

end

