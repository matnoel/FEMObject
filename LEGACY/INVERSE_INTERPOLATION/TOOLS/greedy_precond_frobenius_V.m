function [P,err2,xie,res_F,kappa]=greedy_precond_frobenius_V(A_Xi,mp,P,V,in,constrain)


if nargin<6
    constrain=1;
end

xie=zeros(mp,1);
err2=zeros(mp,1);
kappa={};


res_F(:,1) = in.residual;


for k=1:mp
    disp(k)
    [err2(k),xie(k)] = max(in.residual);
    [P,~,~,in] = frobenius_projection_V_new_point(P,A_Xi,V,xie(k),in,0,constrain);
    res_F = [res_F , in.residual] ;
    if nargout==5
        kappa{k}=compute_condition_number(A_Xi,P);
    end
    
end



