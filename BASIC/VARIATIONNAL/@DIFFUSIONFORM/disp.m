function disp(a)
% function disp(a)

if isempty(a.c1) && isempty(a.c2)
    s = [inputname(1) '(v,u) = int( tau grad(v).grad(u) )'];
elseif isempty(a.c2)
    s = [inputname(1) '(v,u) = int( tau grad(v).c.grad(u) )'];
else
    s = [inputname(1) '(v,u) = int( tau c1.grad(v) c2.grad(u) )'];
end

disp(s)
disp(' ')

