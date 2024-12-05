function [u,result] = solve_onebyone(S,A,b)

if israndom(b)
PC = getPC(b);
else
PC = getPC(A);    
end

u = PCTPMATRIXSUM(PC);
bu = b;
mmmax = getparam(S,'nbfoncmax');
pfixmax=getparam(S,'pfixmax');
pfixtol=getparam(S,'pfixtol');
display_ = getparam(S,'display');
update = getparam(S,'update');

%pfixstagn=getparam(S,'pfixstagn');
tol = getparam(S,'tol');


if display_
   fprintf('       Tensor prod resolution. size = %d, \n',size(A,1)) 
end

    stol = cell(1,mmmax);


errmm = zeros(1,mmmax);
isconv=0;
for mm=1:mmmax
   
l0 = normalize(normalizephi(rands(size(b,1),1,PC)));
l = l0;      
for k=1:pfixmax
for i=1:getnbdim(PC)  
   l{i} =1; 
   mat = simplify(expectnodimmtimes(i,l',A,l));   
   sec = simplify(expectnodimmtimes(i,bu',l));
   if numel(sec)==1 && numel(mat)==1
     l{i} = mat\sec;
   elseif numel(mat)==1
     l{i} = mat\sec;  
   else
       if numel(sec)==1
       sec = sec*mean(PC,i);
       end
       mat = getmatrix(multimtimes(mat',getmasseuni(PC,i)));
       l{i} = mat\sec;
   end
   
   normphi = norm(l{i});
   l{i} = l{i}/normphi;
end

l{0}=1;
l{0} = expect(l,A,l)\expect(bu,l);

err = norm(l-l0)/norm(l);

l0=l;
if display_
fprintf('       iteration %d : erreur = %.3e\n',k,err)
end
    if err<pfixtol
    break
    end
end


stol{mm}=l;

if update && mm>1
AA=zeros(mm,mm);
bb=zeros(mm,1);
for ii=1:mm
    bb(ii)=expect(stol{ii}',b);
    for jj=1:mm
        AA(ii,jj)=expect(stol{ii}',A,stol{jj});
    end
end
ff = AA\bb;
for ii=1:mm
   stol{ii}= ff(ii)*stol{ii};
end

u = PCTPMATRIXSUM(PC);
bu = b;
for ii=1:mm
u = u + stol{ii} ;
bu = bu - A* stol{ii};    
end

else
u = u +l ;
bu = bu - A*l;
end
errmm(mm) = norm(bu)/norm(b);


if errmm(mm)<tol
    if display_ 
    fprintf('       tensor prod resolution : convergence avec nbfun %d : erreur %.3e \n',getm(u),errmm(mm));
    end
    isconv=1;
    break
else
    if display_
    fprintf('       nbfun %d : erreur = %.3e\n',getm(u),errmm(mm));
    end
end
end

if ~isconv
    fprintf('       tensor prod resolution : non-convergence avec nbfun %d : erreur %.3e \n',getm(u),errmm(mm));   
end

result.error = errmm(1:getm(u));
