function ls=setenrich(ls,choix,l)
if nargin==1
    choix=1;
end
if nargin<3
  l = 1:ls.n;  
end
for i=l
ls.LS{i}=setenrich(ls.LS{i},choix);
end