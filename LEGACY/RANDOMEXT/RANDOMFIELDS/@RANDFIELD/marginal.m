function ma=marginal(RF,x)

if nargin==1
ma = RF.marginal;
elseif nargin==2
ma=eval(RF.marginal,x);
end
