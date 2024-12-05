function [P,err2,xie,res]=greedy_update_precond_frobenius_V(A_Xi,mp,P,V,in)


xie=zeros(mp,1);
err2=zeros(mp,1);

for k=1:mp
    disp(k)
    [err2(k),xie(k)] = max(in.residual);
    [P,~,~,in] = frobenius_projection_V_new_point(P,A_Xi,V,xie(k),in,0);
    
end





