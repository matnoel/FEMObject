function [u,result] = separated_decomposition_onebyone(S,b)


if getparam(S,'cyclic')==0
[u,result] = separated_decomposition_onebyone_nocyclic(S,b)    ;  
return
end


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
for k=1:pfixmax
for i=1:getnbdim(PC)  
   l{i} =1;    
   l{i} = simplify(expectnodimmtimes(i,l',l))\simplify(expectnodimmtimes(i,bu',l));
   normphi = norm(l{i});
   l{i} = l{i}/normphi;
end
l{0}=1;
l{0} = expect(l',l)\expect(bu,l);

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
