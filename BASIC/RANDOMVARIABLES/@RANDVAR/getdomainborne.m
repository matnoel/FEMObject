function domainborne = getdomainborne(rv)

domainborne =  getdomain(rv);

if domainborne(1)==-Inf;
    domainborne(1)=min(mean(rv)-5*std(rv));
end

if domainborne(2)==Inf;
    domainborne(2)=max(mean(rv)+5*std(rv));
end

