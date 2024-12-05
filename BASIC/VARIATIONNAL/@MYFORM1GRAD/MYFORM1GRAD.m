function l = MYFORM1GRAD(approx)
% function a =  MYFORM1GRAD(approx)
% a(v,w,u1,u2) = int ( grad(v).grad(w) u1.u2  + grad(v).grad(u1)u2w + grad(v).grad(u2)u1w )

l=struct();
if nargin==0
    l.approx = 0;
else
    l.approx = approx;
end
l = class(l,'MYFORM1GRAD',VARFORM(4));


