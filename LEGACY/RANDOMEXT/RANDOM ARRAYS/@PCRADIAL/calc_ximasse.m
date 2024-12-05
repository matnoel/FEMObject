function R = calc_ximasse(R)
rep = [];
for k=1:R.m
if isempty(R.Lmasse)
rep = [rep,k];    
end
end

if ~isempty(rep)
R.Lmasse(rep) = calc_ximasse(R.POLYCHAOS,R.L(rep)); 
end

