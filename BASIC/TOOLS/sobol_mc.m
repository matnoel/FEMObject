function [Si,V,E] = sobol_mc(X,fun,N,ind)
% function [S,V,E] = sobol_mc(X,fun,N,ind)
% N: nombre de samples
% fun: function handle
% X: RANDVARS
% ind: indices de sobol souhait�s

if nargin<=3
    ind = 1:getM(X);
end

xs1 = random(X,N,1);
xs2 = random(X,N,1);
xs1=[xs1{:}];
xs2=[xs2{:}];

fs1 = zeros(1,N);
fs2 = zeros(1,N);
try
    fs1=fun(xs1);
catch
    for k=1:N
        fs1(k)=fun(xs1(k,:));
    end
end

E = sum(fs1)/N;
V = sum(fs1.^2)/N - E^2;
Vi = zeros(1,length(ind));
for itemp=1:length(ind)
    i = ind(itemp);
    xst = xs2;
    xst(:,i)=xs1(:,i);
    try
        fs2=fun(xst);        
    catch
        for k=1:N
            fs2(k)=fun(xst(l,:));    
        end
    end
    Vi(itemp) = sum(fs1.*fs2)/N - E^2;
end

Si = Vi/V;
