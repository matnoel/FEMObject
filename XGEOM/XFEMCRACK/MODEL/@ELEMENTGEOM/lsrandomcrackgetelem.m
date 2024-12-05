function [elems,reps] = lsrandomcrackgetelem(elem,ls,choix,node)

switch choix
    case 'cut'
if getnbtip(ls)==0
temp = lsispossibly(elem,getlssupport(ls),@lsiscut,node);
else
 lstemp = {getlssupport(ls)};
 for k=1:getnbtip(ls)  
   lstemp=[lstemp,{getlstip(ls,k)}];
 end   
temp = lsispossiblycutnotbicut(elem,lstemp,node);
end
reps = find(temp);
elems = getelem(elem,find(temp));
elems = setlsnumber(elems,getnumber(ls));  
    case 'bicut'
elems = cell(1,getnbtip(ls));
for k=1:getnbtip(ls)        
reps{k} = find(lsispossibly(elem,{getlssupport(ls),getlstip(ls,k)},@lsisbicut,node)) ;  
elems{k} = getelem(elem,reps{k});
end
if getnbtip(ls)==1
elems=elems{1};
reps=reps{1};
end

    case 'in'
[e1,rep1] = lsrandomcrackgetelem(elem,ls,'cut',node);
[e2,rep2] = lsrandomcrackgetelem(elem,ls,'bicut',node);
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

