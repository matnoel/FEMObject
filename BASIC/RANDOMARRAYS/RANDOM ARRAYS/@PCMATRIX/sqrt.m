function x=sqrt(y)

y=calc_masse(y);
x0 = one(getPC(y));
conv=0;
for k=1:100
   x = x0 - 1/2*solve(x0,x0*x0-y);
   %err=norm(x-x0)/norm(x0);
   err = norm(x*x-y)/norm(y);
   %fprintf('sqrt : stagnation error = %f\n',err)
   if err<1e-10
       conv=1;
      break 
   end
   x0 = x;
end

if conv==0
    fprintf('sqrt non converge : erreur %f ',err)
end