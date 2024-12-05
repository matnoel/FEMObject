function e = conditional_expect_samples(fun,i,xi,x,w)
% function e = conditional_expect_samples(fun,i,xi,x,w)

X = 

y = sobol_fun(x,i,xi)

c = [1,2,5,10,20,50,100,500];

if nargin==1
    for i=1:size(x,2)
        x(:,i) = (abs(4*x(:,i)-2)+c(i))/(1+c(i));
    end
    y = prod(x,2);
else
    X = zeros(size(x,1),8);
    j = setdiff(1:8,i);
    X(:,j)=x;
    X(:,i)=xi;
    y = sobol_fun(X);
end
