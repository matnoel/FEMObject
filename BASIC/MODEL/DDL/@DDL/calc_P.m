function P = calc_P(ddl,sys)
nb = length(ddl.ddlgroup);
n = zeros(1,nb);
m = zeros(1,nb);
repn = cell(1,nb);
repm = cell(1,nb);
for i=1:nb
    Pblock{i} = calc_P(ddl.ddlgroup{i},sys);
    n(i) = size(Pblock{i},1);
    m(i) = size(Pblock{i},2);
    repn{i} = sum(n(1:i-1))+[1:n(i)];
    repm{i} = sum(m(1:i-1))+[1:m(i)];
end

P = zerosND([sum(n),sum(m),sizeND(sys)]);

for i=1:nb
    P(repn{i},repm{i}) = Pblock{i};
end