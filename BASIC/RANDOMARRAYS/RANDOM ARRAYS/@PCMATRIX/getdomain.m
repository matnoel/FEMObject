function domain = getdomain(a)

r = double(random(a,1e3));

domain = [min(r,[],2),max(r,[],2)];
domainborne = domain ;
