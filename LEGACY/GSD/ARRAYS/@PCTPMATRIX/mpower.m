function x = mpower(x,i)

if i~=int8(i)
    error('rentrer une puissance entiere')
end

x0 = x;
for l=2:i
x = mtimes(x0,x);
end
