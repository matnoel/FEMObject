function [P,residual]=greedy_inverse_interp_frob(A,Y,m,V)
% function [P,residual,out]=greedy_inverse_interp_frob(A,P,m,V,opts)
%
% Greedy selection of interpolation poits
% A : operator whose inverse is approached
% P : initial precond (possibly P=[])
% m : number of (new) interpolation points
% opts : options

% if nargin<4 || isempty(V)
%     V=eye(A.sz(1,:));
% end
% 
% if nargin<2 || isempty(P)
%     
% end

residual=cell(m+1,1);
xie=zeros(m,1);

opts.display=0;
[P,residual{1},out] = inverse_interp_frob(A,Y,V,opts);

PivotElementM_1 = out.PivotElementM_1;
n_evalM = length(PivotElementM_1);
PivotElementS_1 = out.PivotElementS_1;
n_evalS = length(PivotElementS_1);
PAV_M = out.PAV_M;
M_Phi = out.M_Phi;
S_Phi = out.S_Phi;
M_xi = out.M_xi;
S_xi = out.S_xi;

clear out

for k=1:m
    fprintf('%d ',k)
    [~,xie(k)] = max(residual{k});
    
    Y_new=implicit_inverse_lu( eval_param_operator(A,xie(k)) );
    Y = [Y;{Y_new}];
    
    % update M
    for i=1:n_evalM
        PAV_new = Y_new * ( eval_param_operator(A,PivotElementM_1(i))*V );
        PAV_new = PAV_new(:);
        
        b = ( PAV_new'*PAV_M{i})';
        a = PAV_new'*PAV_new;
        M_xi{i} = [M_xi{i},b ; b' , a];
        
        PAV_M{i} = [PAV_M{i} , PAV_new];
    end
    
    % update S
    parfor i=1:n_evalS
        PAV_new = Y_new * ( eval_param_operator(A,PivotElementS_1(i))*V );
        S_xi{i} = [S_xi{i} ; trace( V'*PAV_new )];
    end
    
    % projection    
    [Phi,residual{k+1}]=MS_solve(M_xi,M_Phi,S_xi,S_Phi);
    residual{k+1}=residual{k+1} + trace(V'*V);
    residual{k+1}(residual{k+1}<0)=0;
    
    
end

PT = TSPACE_OPERATORS( [{Y};{param_operator( Phi )}] );
PC = CANONICAL_CORE( ones(PT.dim(1),1),2 );
P = LRTENSOR(PC,PT);


