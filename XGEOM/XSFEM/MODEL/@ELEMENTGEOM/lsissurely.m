function issurely = lsissurely(elem,ls,fun,node)
if ~isa(ls,'cell')
    ls={ls};
end
lsxnode=cell(1,length(ls));
for i=1:length(ls)
if nargin==4
lsxnode{i} = getconnecvalue(elem,ls{i},node);
else
lsxnode{i}=ls{i};
end
end


if length(ls)==1

issurely = fun(elem,lsxnode{1}(:,:,1));
for i=2:size(lsxnode{1},3)
issurely = issurely & fun(elem,lsxnode{1}(:,:,i));
end

elseif length(ls)==2
    
issurely = fun(elem,lsxnode{1}(:,:,1),lsxnode{2}(:,:,1));
for i=2:max(size(lsxnode{1},3),size(lsxnode{2},3))
    if size(lsxnode{1},3)==1
        i1=1;
    else
        i1=i;
    end
    if size(lsxnode{2},3)==1
        i2=1;
    else
        i2=i;
    end
issurely = issurely & fun(elem,lsxnode{1}(:,:,i1),lsxnode{2}(:,:,i2));
end

else
    error('pas programme')
 
end





