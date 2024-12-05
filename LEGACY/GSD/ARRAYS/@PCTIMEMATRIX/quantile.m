function u = quantile(u,q,n)
% function v = quantile(u,q,n)
% u : PCTIMEMATRIX
% q : quantiles
% n : nombre de tirages pour le calcul des quantiles
% v : TIMEMATRIX

if nargin==2
    n=1000;
end


if isa(u.value,'cell') 
   u = cell2mat(u);
end

   
   if ~all(size(u)==1)
       error('pas programme pour plus d''une composante')
   end
   u.value = double(random(u.value,n));
   u.value = quantile(u.value',q);
   
   u.s=[length(q),1];
   u = TIMEMATRIX(u.value,u.TIMEMODEL,u.s);

   
