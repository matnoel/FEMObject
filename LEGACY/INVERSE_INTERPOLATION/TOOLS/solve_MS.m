function [Phi,residual]=solve_MS(M,S)



n_eval=M.n;

if M.n>100
    Phi = cell(n_eval,1);
    residual = zeros(n_eval,1);
    
    
    parfor i=1:n_eval
        MM = eval(M,i);
        SS = eval(S,i);
        lambda = MM\SS;
        Phi{i} = lambda;
        residual(i)= lambda'*MM*lambda - 2*lambda'*SS;
    end
    
    
else
    Phi = cellfun(@(ii) eval(M,ii)\eval(S,ii)...
     , num2cell(1:M.n) , 'UniformOutput',0);
    residual = cellfun(@(ii,lambda) lambda'*eval(M,ii)*lambda - 2*lambda'*eval(S,ii)  ...
                 , num2cell(1:M.n),Phi);
end



