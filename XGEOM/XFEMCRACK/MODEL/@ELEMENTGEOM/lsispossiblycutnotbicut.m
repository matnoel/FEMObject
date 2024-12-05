function ispossibly = lsispossiblycutnotbicut(elem,ls,node)
if ~isa(ls,'cell')
    error('renter les levelsets sous forme de cellule')
end
lsxnode=cell(1,length(ls));
for i=1:length(ls)
if nargin==3
lsxnode{i} = getconnecvalue(elem,ls{i},node);
else
lsxnode{i} = ls{i};
end
end

istemp = lsiscutnotbicut(elem,lsxnode{1}(:,:,1),lsxnode{2}(:,:,1));
for k=3:length(lsxnode)
istemp = istemp & lsiscutnotbicut(elem,lsxnode{1}(:,:,1),lsxnode{k}(:,:,1));    
end
ispossibly = istemp;
imax = 1;
for k=1:length(lsxnode)
imax=max(imax,size(lsxnode{k},3));    
end

for i=2:imax
    for k=1:length(lsxnode)    
    if size(lsxnode{k},3)==1
        ik{k}=1;
    else
        ik{k}=i;
    end
    end
    
istemp = lsiscutnotbicut(elem,lsxnode{1}(:,:,ik{1}),lsxnode{2}(:,:,ik{2}));
for k=3:length(lsxnode)
istemp = istemp & lsiscutnotbicut(elem,lsxnode{1}(:,:,ik{1}),lsxnode{k}(:,:,ik{k}));
end
ispossibly = ispossibly | istemp;
end    


