function [a,xi] = randomblock(x,n,m)

if n<=m
    [a,xi] = random(x,n,1);
else    
    [a,xi] = random(x,m,1);    
    j = n-m;
    while j>0
        pourcentage(n-j,n);
        [atemp,xitemp] = random(x,min(j,m),1);  
        a = [a;atemp];
        xi = [xi;xitemp];
        j = j - min(j,m);
    end
end


