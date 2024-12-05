function l = MYFORM1(approx)
% function a =  MYFORM1()
% a(v,u) = int u^2 grad(u) grad(v)

l=struct();
if nargin==0
    l.approx = 0;
else
    l.approx = approx;
end

l = class(l,'MYFORM1',VARFORM(2));


