function P = shepard_precond(P,Xi,Xie,s)



nXie=size(Xie,1);

LBD_S=zeros(P.r,size(Xi,1));
for k=1:size(Xi,1)
    
    dk = ( sum( (Xie - repmat( Xi(k,:) ,[nXie,1]) ).^2 ,2)).^(0.5)  +eps;
    lambda = dk.^(-s)/sum(dk.^(-s));
    
    LBD_S(:,k)= lambda;
    
end


P.Phi=LBD_S;
P.n=size(Xi,1);