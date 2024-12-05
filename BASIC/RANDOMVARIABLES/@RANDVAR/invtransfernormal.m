function ax = invtransfernormal(a,xix)
% function ax = invtransfernormal(a,xix)
% a : RANDVAR    
% ax=F_a^-1(Phi(xix))
xi=RVNORMAL;
ax=transfer(xi,a,xix);
