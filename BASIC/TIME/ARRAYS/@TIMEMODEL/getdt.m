function dt = getdt(T,i)
% function dt = getdt(T,i)

if nargin==1
    dt=T.dt;
else
   dt = T.dt(i); 
end

