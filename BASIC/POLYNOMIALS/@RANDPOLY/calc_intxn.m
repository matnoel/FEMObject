function intxn=calc_intxn(h,liste)

param = getparam(h);
if isfield(param,'intxn') && (length(param.intxn)-1 >= max(liste)) 
    intxn = param.intxn(liste+1);
else
    
intxn=zeros(1,length(liste));
param = getparam(h);
    r=param.r;
    d = getdomain(r);
  for i=1:length(liste)
 n=liste(i);
%  intxn(i)=quad(@(x)evaldensity(h,x).*x.^n,d(1),d(2),eps,0); 
  intxn(i)=integral(@(x)evaldensity(h,x).*x.^n,d(1),d(2),'AbsTol',eps);   end
end

