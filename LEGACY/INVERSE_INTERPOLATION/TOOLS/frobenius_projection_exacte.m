function [P,M,S,out] = frobenius_projection_exacte(P,A,display,constrain)

if nargin<3
    display=1;
end

if nargin<4
    constrain=1;
end

%% Compute M_ij(\xi) = trace( A(\xi)^t P_i^t P_j A(\xi) )
if display
disp('Compute M_ij(\xi) = trace( A(\xi)^t P_i^t P_j A(\xi) )')
end
tol = 1e-15;

% Compute the Phi function of the product A*A
[I,J]=ind2sub([A.r,A.r],1:(A.r*A.r));
Phi = A.Phi(I,:).*A.Phi(J,:);

% Cross approx of Phi : new function Phi, and evaluation point PivotElement
[PivotValue,PivotElement,Rows,Cols,ifail] = CompleteACA(Phi',tol);
PivotElement=PivotElement(:,1);
n_eval = length(PivotElement);
Phi = (Cols/Cols(PivotElement,:))';



MA=cell(n_eval,1);


parfor i=1:n_eval
    
    if display
    fprintf('\n%d/%d ', i,n_eval)
    end
    
    A_eval = eval_sparse(A,PivotElement(i)) ;
    tmp=cellfun( @(L,U) U\(L\full(A_eval)) ,P.L,P.U,'UniformOutput',0);
    tmp=cellfun( @(x) x(:),tmp,'UniformOutput',0);
    tmp=[tmp{:}];
    
    MA{i} = tmp'*tmp;
end

M = affine_matrix(MA,Phi);


%% Compute S_i(\xi) = trace( P_i A(\xi) )
if display
disp('Compute S_i(\xi) = trace( P_i A(\xi) )')
end


Phi = A.Phi;
% Cross approx of Phi : new function Phi, and evaluation point PivotElement
[PivotValue,PivotElement,Rows,Cols,ifail] = CompleteACA(Phi',tol);
PivotElement=PivotElement(:,1);
n_eval = length(PivotElement);
Phi = (Cols/Cols(PivotElement,:))';

SA=cell(n_eval,1);

parfor i=1:n_eval
    
    if display
    fprintf('\n%d/%d ', i,n_eval)
    end
    
    A_eval = eval_sparse(A,PivotElement(i)) ;
    tmp=cellfun( @(L,U) U\(L\full(A_eval)) ,P.L,P.U,'UniformOutput',0);
    
    SA{i}=cellfun(@trace,tmp);
    
    
    
end


S = affine_matrix(SA,Phi);



%% Compute M\S
if display
disp('Compute M\S')
end
if constrain
    [Phi,residual]=solve_MS_const(M,S);
else
    [Phi,residual]=solve_MS(M,S);
end

residual=residual + A.s(1);
residual(residual<0)=0;
out.residual=sqrt(residual);

P.Phi = [Phi{:}];
P.n = size(P.Phi,2);



