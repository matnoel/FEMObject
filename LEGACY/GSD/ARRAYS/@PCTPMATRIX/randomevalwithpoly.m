function a = randomevalwithpoly(x,Hs)

phi = ones(size(Hs{1},1),1);
isr = isranddim(x);
for i=1:getnbgroups(x)
    if isr(i)
    phi = phi.*(Hs{i}*x.phi{i});
    else
    phi = phi.*x.phi{i};    
    end
end

if numel(phi)==1 || numel(x)==1
    a = x.phi0*phi;    
else
    a =  MULTIMATRIX(x.phi0(:)*phi(:)',size(x),[numel(phi),1]);
end
    
