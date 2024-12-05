function [Phi,residual]=MS_solve(M_xi,M_Phi,S_xi,S_Phi)


s=size(M_xi{1},1);
M_xi = cellfun(@(x) x(:), M_xi,'UniformOutput',0);
M_xi = [M_xi{:}];
S_xi = [S_xi{:}];

n_eval = size(M_Phi,1);
Phi = cell(n_eval,1);
residual = zeros(n_eval,1);


parfor i=1:n_eval
    
    MM = M_xi*(M_Phi(i,:)');
    MM = reshape(MM,s,s);
    SS = S_xi*(S_Phi(i,:)');
    lambda = MM\SS;
    
    Phi{i} = lambda;
    residual(i)= lambda'*MM*lambda - 2*lambda'*SS;
end

Phi = [Phi{:}];
