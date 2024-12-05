function x = heaviside_filtered(x,l,a,c)
% function x = heaviside_filtered(x,l,a,c)
% H(x) = 1 - signe(x)*exp(-c*|x/l|^a)
% l = (max(x)-min(x))/100 par d�faut
% a=1 et c=5 par d�faut

if l==0
    x = double(x>0);
else
    if nargin<=2
        a=1;
    end
    if nargin<=3
        c = 5;
    end
    if nargin<=1
        l = (max(max(x))-min(min(x)))/100;
    end

%x = (x>0).*(1-1/2*exp(-c*abs(x/l).^a))+(x<=0).*(1/2*exp(-c*abs(x/l).^a));
    x = (x<2*l).*1/2.*(1+tanh(3*x/l)) + (x>2*l);
end


