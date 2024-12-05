function hm=mean(h,liste)

hm = sparse(length(liste),1);
hm(find(liste==0))=1;