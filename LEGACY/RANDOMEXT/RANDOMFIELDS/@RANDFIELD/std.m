function s=std(RF,x)

s=std(eval(RF.marginal,x));

if numel(s)==1
s=repmat(s,size(x,1),1);
end