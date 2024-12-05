function [P,xie,res]=greedy_delta(A,Xi,U,m,H,constrain,b)



if nargin<6
    constrain=1;
end

xie       = zeros(m+1,1);
borne     = cell(m+1,1);
delta     = cell(m+1,1);
error     = cell(m+1,1);
borne_sup = zeros(m+1,1);


P = implicit_affine_matrix( {speye(A.s)} ,ones(1,A.n));
[delta{1},borne{1},error{1}] = compute_delta(A,P,U,b);

[borne_sup(1),xie(1)] = max(borne{1});
% res_F(:,1) = in.residual;


for k=1:m
    disp(k)
    
    
    if k==1
        P = implicit_affine_matrix( {eval_sparse(A,xie(1))} ,1);
        [P,~,~,in] = frobenius_projection_V(P,A,H,0,constrain);
    else
        [P,~,~,in] = frobenius_projection_V_new_point(P,A,H,xie(k),in,0,constrain);
    end
    
    [delta{k+1},borne{k+1},error{k+1}] = compute_delta(A,P,U,b);
    
    [borne_sup(k+1),xie(k+1)] = max(borne{k+1});
    
    
    
    
    
    
end


res.borne_sup=borne_sup;
res.borne=borne;
res.delta=delta;
res.stats=in;
res.error=error;