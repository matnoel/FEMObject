function [P,M,S,out] = frobenius_projection_V(P,A,V,display,constrain)

if nargin<4
    display=1;
end

if nargin<5
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% salut=1;
% 
% PAV = cellfun( @(ii) eval_sparse(A,ii)*V , num2cell(PivotElement) ,'UniformOutput',0) ;
% PAV = cellfun( @(AV) cellfun(@(L,U) U\(L\AV),P.L,P.U, 'UniformOutput',0) , PAV ,'UniformOutput',0) ;
% PAV = cellfun( @(PAV) cellfun(@(x) x(:) ,PAV, 'UniformOutput',0) , PAV ,'UniformOutput',0) ;
% PAV = cellfun( @(PAV) cellfun(@(x) x(:) ,PAV, 'UniformOutput',0) , PAV ,'UniformOutput',0) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


MA=cell(n_eval,1);
nargout4=(nargout == 4);
if nargout4
    out.PAV_M=cell(n_eval,1);
end

parfor i=1:n_eval
    
    if display
    fprintf('\n%d/%d ', i,n_eval)
    end
    
    Axi = eval_sparse(A,PivotElement(i)) ;
    
    AV=Axi*V;
    
    PAV = cellfun(@(L,U) U\(L\AV),P.L,P.U, 'UniformOutput',0);
    
    if nargout4
        PAV_M{i}=PAV;
    end
    
    PAV = cellfun(@(x) x(:),PAV, 'UniformOutput',0);
    PAV = [PAV{:}];
    MA{i} = PAV'*PAV;
    
    
end

M = affine_matrix(MA,Phi);


if nargout == 4
    out.PAV_M=PAV_M;
    out.Phi_M = Phi;
    out.PivotElement_M = PivotElement;
    out.M = M;
end


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
    Axi = eval_sparse(A,PivotElement(i)) ;
    AV=Axi*V;
    VPAV = cellfun(@(L,U) V'*( U\(L\AV) ),P.L,P.U, 'UniformOutput',0);
    SA{i} = cellfun(@(x)trace(x),VPAV, 'UniformOutput',1);
    
end


S = affine_matrix(SA,Phi);

if nargout4
    out.Phi_S = Phi;
    out.PivotElement_S = PivotElement;
    out.S = S;
end


%% Compute M\S
if display
disp('Compute M\S')
end
if constrain
    [Phi,residual]=solve_MS_const(M,S);
else
    [Phi,residual]=solve_MS(M,S);
end

residual=residual + trace(V'*V);
residual(residual<0)=0;
out.residual=sqrt(residual);

P.Phi = [Phi{:}];
P.n = size(P.Phi,2);




