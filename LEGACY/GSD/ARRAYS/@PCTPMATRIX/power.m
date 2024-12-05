function x = power(x,i)

if ~isinteger(i)
    error('rentrer une puissance entiere')
end

x0 = x;
for l=2:i
x = times(x0,x);
end
