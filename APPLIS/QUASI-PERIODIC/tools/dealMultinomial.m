function repartition = dealMultinomial(probability,dealSize)
% repartition = dealMultinomial(probability,dealSize)

deal = mnrnd(1,probability,dealSize) ;

assert(~isnan(deal(1)),'Sum of probabilities must be one')

dealNb = size(deal,2) ;
repartition = cell(dealNb,1) ;
for v = 1:dealNb
    repartition{v} = find(deal(:,v)) ;
end

end