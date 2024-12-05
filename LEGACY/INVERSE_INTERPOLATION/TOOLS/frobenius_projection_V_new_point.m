function [P,M,S,in] = frobenius_projection_V_new_point(P,A,V,xi,in,display,constrain)

if nargin<6
    display=1;
end
if nargin<7
    constrain=1;
end

[P.L{end+1,1},P.U{end+1,1}]=lu(eval_sparse(A,xi));
P.r=P.r+1;
P.n=P.n+1;

%% Compute M_ij(\xi) = trace( A(\xi)^t P_i^t P_j A(\xi) )
if display
    disp('Compute M_ij(\xi) = trace( A(\xi)^t P_i^t P_j A(\xi) )')
end


% Load the Phi function of the product A*A
Phi = in.Phi_M;
% Load PivotElement
PivotElement=in.PivotElement_M;
n_eval = length(PivotElement);


parfor i=1:n_eval
    PAV{i} = eval_sparse(A,PivotElement(i))*V;
    PAV{i} = P.U{end}\( P.L{end}\PAV{i} );
end

% PAV = cellfun( @(ii) eval_sparse(A,ii)*V , num2cell(PivotElement) ,'UniformOutput',0);
% PAV = cellfun( @(AV) P.U{end}\( P.L{end}\AV ) , PAV ,'UniformOutput',0) ;
in.PAV_M = cellfun(@(a,b) [a;{b}],in.PAV_M,PAV,'UniformOutput',0);

PAV = cellfun( @(PAV) cellfun(@(x) x(:) ,PAV, 'UniformOutput',0) , in.PAV_M ,'UniformOutput',0) ;
PAV = cellfun( @(PAV) [PAV{:}] , PAV ,'UniformOutput',0) ;
PAV = cellfun( @(PAV) PAV'*PAV , PAV ,'UniformOutput',0) ;

M = affine_matrix(PAV,Phi);
in.M=M;

%% Compute S_i(\xi) = trace( P_i A(\xi) )
if display
    disp('Compute S_i(\xi) = trace( P_i A(\xi) )')
end

% Load evaluation point PivotElement
PivotElement=in.PivotElement_S;
n_eval = length(PivotElement);
Phi = in.Phi_S;

salut=1;

PAV = cellfun( @(ii) eval_sparse(A,ii)*V , num2cell(PivotElement) ,'UniformOutput',0);
PAV = cellfun( @(x) V'*( P.U{end}\( P.L{end}\x ) ) , PAV ,'UniformOutput',0) ;
PAV = cellfun( @(x) trace(x) , PAV ,'UniformOutput',0) ;

SA=cellfun(@(a,b)[a;b], in.S.A , PAV ,'UniformOutput',0);

S = affine_matrix(SA,Phi);
in.S=S;


%% Compute M\S
if constrain
    [Phi,residual]=solve_MS_const(M,S);
else
    [Phi,residual]=solve_MS(M,S);
end

residual=residual + trace(V'*V);
residual(residual<0)=0;
in.residual=sqrt(residual);

P.Phi = [Phi{:}];
P.n = size(P.Phi,2);


% %%% MIN sans contrainte de positivite
% Phi = cellfun(@(ii) eval(M,ii) \ eval(S,ii) , num2cell(1:A.n) , 'UniformOutput',0);
% P.Phi = [Phi{:}];




