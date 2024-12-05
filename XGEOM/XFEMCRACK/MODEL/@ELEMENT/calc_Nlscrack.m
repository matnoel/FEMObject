function Nls=calc_Nlscrack(elem,xnode,xgauss,ls,choix)

lssupport = getlssupport(ls);
ls1 = lssupport(getconnec(elem)');

connecenrich = getparam(elem,'connecenrich');
connecnbddl = getparam(elem,'connecnbddl');

if nargin==5
switch choix
    case 'in'
      fact=-1;  
    case 'out'
      fact=1;
end
else
ls1gauss = getN(elem,xgauss)*ls1;
fact = (ls1gauss>=0)-(ls1gauss<0);
end


N = calc_N(elem,xnode,xgauss);
n=size(N,1);
Nls = zerosND(n,getnbddl(elem),getnbelem(elem),size(N,4)); 

for e=1:getnbelem(elem)
for i=1:3
repddli = n*(i-1)+[1:n];
Ni = N(:,repddli,e); 
if connecenrich(e,i)
Nlsi = [Ni , fact.*Ni];
else
Nlsi=Ni;
end
repddlienrich = sum(connecnbddl(e,1:i-1))+[1:connecnbddl(e,i)];
try
Nls(:,repddlienrich,e,:)=Nlsi;
catch
   keyboard 
end
end
end

