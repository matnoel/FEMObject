function PC=calc_P(PC)

if ~isempty(PC.indices)
    PC.P = size(PC.indices,1)-1;
elseif PC.typebase==2
    PC.P = prod(PC.n)-1;
elseif length(unique(PC.n))==1
    temp = unique(PC.n)-1;
    PC.P = factorial(temp+PC.M)/factorial(temp)/factorial(PC.M)-1;
else
    PC.P=[];
    warning('pas de calcul de P')
end
