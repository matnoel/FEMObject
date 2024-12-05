function intxn=calc_intxn(h,liste)

intxn=zeros(1,length(liste));

for i=1:length(liste)    
    n=liste(i);
    param=get(h,'param');
    a=param.a;
    intxn(i)=gamma(a+n)/gamma(a);

end
