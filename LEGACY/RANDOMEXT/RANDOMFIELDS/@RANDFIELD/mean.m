function mu=mean(RF,x)

mu=mean(eval(RF.marginal,x));

if numel(mu)==1
mu=repmat(mu,size(x,1),1);
end