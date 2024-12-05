function u=reshapeall(u,s,sm)
if nargin==3
   if prod(s)*prod(sm)~=prod(u.s)*prod(u.sm)
    error('To RESHAPE the number of elements must not change.')
   end
   u.s = s;
   u.sm = sm;
   
else
    error('enter  3 arguments')
end

    