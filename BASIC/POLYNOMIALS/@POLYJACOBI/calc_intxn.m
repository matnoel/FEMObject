function intxn=calc_intxn(h,liste)
%
intxn=zeros(1,length(liste));
param = getparam(h);
    a=param.a;
    b=param.b;  
intex=1;
if intex
  for i=1:length(liste)
 n=liste(i);
  intxn(i)=integral(@jacobixn,-1,1,eps,0,n,a,b); 
 end
else
gauss = calc_gausspoints(POLYLEGENDRE(),floor((max(liste)+1)/2)+15);
for i=1:length(liste)    
n=liste(i);
intxn(i)=2*gauss.w * jacobixn(gauss.coord,n,b,a);
end
end

function y=jacobixn(x,n,b,a)

y=(x.^n).*((1+x).^(a-1)).*((1-x).^(b-1))/beta(a,b)/2^(a+b-1);

return