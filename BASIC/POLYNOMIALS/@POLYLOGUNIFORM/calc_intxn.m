function intxn=calc_intxn(h,liste)
%

param = getparam(h);
r = param.r;
A = getparam(r,'A');
B = getparam(r,'B');

intxn=zeros(1,length(liste));
rep = find(liste~=0);
intxn(rep) = (B.^liste(rep)-A.^liste(rep))./liste(rep)/(log(B)-log(A));
rep = find(liste==0);
intxn(rep)=1;

 
 