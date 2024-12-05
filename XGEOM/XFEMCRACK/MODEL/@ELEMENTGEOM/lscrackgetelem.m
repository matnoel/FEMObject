function [elems,reps] = lscrackgetelem(elem,ls,choix,node)


switch choix
    case 'cut'
if getnbtip(ls)==0
temp = lsiscut(elem,getlssupport(ls),node);
else
temp = lsiscutnotbicut(elem,getlssupport(ls),getlstip(ls,1),node);  
for k=2:getnbtip(ls)
    temp = temp & lsiscutnotbicut(elem,getlssupport(ls),getlstip(ls,k),node);  
end
end
reps = find(temp);
elems = getelem(elem,find(temp));
elems = setlsnumber(elems,getnumber(ls));  
    case 'bicut'
elems = cell(1,getnbtip(ls));
for k=1:getnbtip(ls)        
reps{k} = find(lsisbicut(elem,getlssupport(ls),getlstip(ls,k),node)) ;  
elems{k} = getelem(elem,reps{k});
end
if getnbtip(ls)==1
elems=elems{1};
reps=reps{1};
end


    case 'in'
[e1,rep1] = lscrackgetelem(elem,ls,'cut',node);
[e2,rep2] = lscrackgetelem(elem,ls,'bicut',node);
reps = rep1 ; 
if isa(rep2,'cell')
    for k=1:length(rep2)
        reps = union(reps,rep2{k});
    end
else
    reps = union(reps,rep2);
end

reps = setdiff(1:getnbelem(elem),reps);
elems = getelem(elem,reps);

    otherwise
        error('choix pas bon : cut ou bicut')
end
%keyboard


