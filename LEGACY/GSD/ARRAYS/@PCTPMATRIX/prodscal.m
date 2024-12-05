function n = prodscal(x,y)

n = sum(sum(expectnodimtimes([],x,y)));
n = simplify(n);
