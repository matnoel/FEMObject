function PC=calc_metric(PC)
% function metric=calc_metric(PC,PC2)
% PC : POLYCHAOS
% metric(alpha,beta)=E(H_alpha H_beta)

ind = getindices(PC);

PC = calc_metricuni(PC);
mpc = getmetricuni(PC);

k=1;
metric = mpc{k}(ind(:,k)+1,ind(:,k)+1);
for k=2:getM(PC)
   metric = metric.*mpc{k}(ind(:,k)+1,ind(:,k)+1);
end

PC.metric = metric;