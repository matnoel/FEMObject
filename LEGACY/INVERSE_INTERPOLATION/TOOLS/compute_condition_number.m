function [kappa]=compute_condition_number(A_Xi,P)

kappa=zeros(A_Xi.n,1);


parfor k_eval=1:A_Xi.n
%     k_eval
    
    A_eval = eval_sparse(A_Xi,k_eval);
    tmp=cellfun( @(L,U) U\(L\full(A_eval)) ,P.L,P.U,'UniformOutput',0);
    
    % Nearest
    PAk=tmp{1}*P.Phi(1,k_eval) ;
    for k=2:P.r
        PAk=PAk+tmp{k}*P.Phi(k,k_eval) ;
    end
    kappa(k_eval)=my_cond(PAk);
    
    
end






