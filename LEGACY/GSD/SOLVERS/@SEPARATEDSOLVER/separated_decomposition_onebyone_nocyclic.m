function [u,result] = separated_decomposition_onebyone_nocyclic(S,b)


PC = getPC(b);

u = PCTPMATRIXSUM(PC);
bu = b;
mmmax = getparam(S,'nbfoncmax');
pfixmax=getparam(S,'pfixmax');
pfixtol=getparam(S,'pfixtol');
%pfixstagn=getparam(S,'pfixstagn');
tol = getparam(S,'tol');

errmm = zeros(1,mmmax);
for mm=1:mmmax
   
l0 = normalize(normalizephi(rands(size(b,1),1,PC)));
l = l0;     
lplus = l;   

for k=1:pfixmax
    
ltemp=l0;    
    
for kkk=0:getnbdim(PC)  

for i=1:getnbdim(PC)  
    
     l = ltemp;
    
   l{i} =1;    
   lplus{i} = simplify(expectnodimmtimes(i,l',l))\simplify(expectnodimmtimes(i,bu',l));
   normphi = norm(lplus{i});
   lplus{i} = lplus{i}/normphi;
end

l=ltemp;
l{0}=1;
lplus{0} = expect(l',l)\expect(bu,l);

ltemp=lplus;

end
l=lplus;

err = norm(l-l0)/norm(l);

l0=l;

fprintf('iteration %d : erreur = %.3e\n',k,err)
    if err<pfixtol
    break
    end
end
u = u +l ;
bu = bu - l;
errmm(mm) = norm(u-b)/norm(b);


if errmm(mm)<tol
    fprintf('convergence avec nbfun %d : erreur %.3e \n',getm(u),errmm(mm));
    break
else
    fprintf('nbfun %d : erreur = %.3e\n',getm(u),errmm(mm));
end
end


result.error = errmm(1:getm(u));
