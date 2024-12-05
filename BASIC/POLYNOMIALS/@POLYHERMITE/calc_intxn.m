function intxn=calc_intxn(h,liste)
%

intxn=zeros(1,length(liste));

for i=1:length(liste)    
    n=liste(i);
    if mod(n,2)==0
    intxn(i)=prod(n-1:-2:1);
    end
end
