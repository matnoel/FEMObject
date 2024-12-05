function ux = eval_sol(S,u,x)
% function ux = eval_sol(S,u,x)
if S.dim==1
   Lx = polyval(S.L,0:getnbpoints(S.L)-1,x);
   ux= Lx*u;

else
nx = getnbpoints(S.L);
ny = getnbpoints(S.L);
[I,J]=meshgrid(1:nx,1:ny);
   Lx = polyval(S.L,0:nx-1,x(:,1));
   Ly = polyval(S.L,0:ny-1,x(:,2));
   ux= (Lx(:,I(:)).*Ly(:,J(:)))*u;
   
end