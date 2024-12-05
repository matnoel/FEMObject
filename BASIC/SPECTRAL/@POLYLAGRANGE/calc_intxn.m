function intxn=calc_intxn(h,liste)

intxn=zeros(1,length(liste));
d = getdomain(h);
for i=1:length(liste)    
    n=liste(i);
    %if mod(n,2)==0
   intxn(i)=(d(2)^(n+1)-d(1)^(n+1))/(n+1)/(d(2)-d(1)); 
    %end
end

