function H = intersect(H1,H2)
% function H = intersect(H1,H2)

H=H1 ;

[rep,b] = ismember(H1,H2);
a=find(rep);
for j=1:length(a)
    H.h{a(j)}=  intersect(H1.h{a(j)},H2.h{b(j)});
end

