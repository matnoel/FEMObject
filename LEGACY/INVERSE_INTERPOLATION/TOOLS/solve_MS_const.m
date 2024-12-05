function [Phi,residual]=solve_MS_const(M,S)



options.Display = 'off';
options.Algorithm = 'trust-region-reflective';

n_eval = M.n;
if M.n>50
    Phi = cell(n_eval,1);
    residual = zeros(n_eval,1);
    
    
%     isOpen = matlabpool('size') > 0;
%     if ~isOpen
%         for i=1:n_eval
%             MM=eval(M,i);
%             MM=(MM+MM')/2;
%             SS=eval(S,i);
%             lambda = quadprog(MM,-SS,[],[],[],[], zeros(M.s(1),1),[],[],options );
%             Phi{i} = lambda;
%             residual(i)= lambda'*MM*lambda - 2*lambda'*SS;
%         end
%         
%     else
        parfor i=1:n_eval
            MM=eval(M,i);
            MM=(MM+MM')/2;
            SS=eval(S,i);
            lambda = quadprog(MM,-SS,[],[],[],[], zeros(M.s(1),1),[],[],options );
            Phi{i} = lambda;
            residual(i)= lambda'*MM*lambda - 2*lambda'*SS;
        end
%     end
    
else
    Phi = cellfun(@(ii) quadprog((eval(M,ii)+eval(M,ii)')/2,-eval(S,ii),[],[],[],[], zeros(M.s(1),1),[],[],options ) ...
        , num2cell(1:M.n) , 'UniformOutput',0);
    
    residual = cellfun(@(ii,lambda) lambda'*eval(M,ii)*lambda - 2*lambda'*eval(S,ii)  ...
        , num2cell(1:M.n),Phi);
end
