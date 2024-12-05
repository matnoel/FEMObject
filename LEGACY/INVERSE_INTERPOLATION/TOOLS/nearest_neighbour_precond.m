function P = nearest_neighbour_precond(P,Xi,Xie)



nXie=size(Xie,1);
LBD_N=zeros(P.r,size(Xi,1));
for k=1:size(Xi,1)
    
    [~,ind]=min(( sum( (Xie - repmat( Xi(k,:) ,[nXie,1])).^2 ,2)));
    LBD_N(ind,k) = 1;
    
end


P.Phi=LBD_N;
P.n=size(Xi,1);