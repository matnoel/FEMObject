function [Phi,residual]=solve_MS_const_K(M,S,gamma_m,gamma_p,C,K)



options.Display = 'off';
options.Algorithm = 'active-set';

n_eval = M.n;
if M.n>50
    Phi = cell(n_eval,1);
    residual = zeros(n_eval,1);
    
    
    parfor i=1:n_eval
        MM=eval(M,i);
        MM=(MM+MM')/2;
        MM=[MM,-MM;-MM,MM];
        SS=eval(S,i);
        SS=[SS;-SS];
        CONST=[ gamma_m , -gamma_p ; K*gamma_m-C , -K*gamma_p-C ];
        
        lambda = quadprog(MM,-SS,-CONST,[0;0],[],[], zeros(size(MM,1),1),[],[],options );
        Phi{i} = [eye(size(MM,1)/2),-eye(size(MM,1)/2)]*lambda;
        residual(i)= lambda'*MM*lambda - 2*lambda'*SS;
    end
    
else
    Phi = cellfun(@(ii) quadprog((eval(M,ii)+eval(M,ii)')/2,-eval(S,ii),[],[],[],[], zeros(M.s(1),1),[],[],options ) ...
        , num2cell(1:M.n) , 'UniformOutput',0);
    
    residual = cellfun(@(ii,lambda) lambda'*eval(M,ii)*lambda - 2*lambda'*eval(S,ii)  ...
        , num2cell(1:M.n),Phi);
end
