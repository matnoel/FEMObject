function w=multisum(u,k)

if nargin==1
if isa(u.value,'cell')
w = u.value{1};
for k=2:numel(u.value)
w = w + u.value{k};    
end
else
    w=reshape(sum(u.value,2),u.s);
end
else
if isa(u.value,'double')
w = switchmulti(sum(switchmulti(u),k));
if prod(sizem(u))==1
w = reshape(w.value,w.s);  
end
else
w=u;
if length(u.sm)>2
    error('pas prevu')
end

if k==1
    w.sm = [1,size(u.value,2)];
    w.value = cell(1,size(u.value,2));
    for kk=1:size(u.value,2)
    w.value{kk} = u.value{1,kk};    
        for ll=1:size(u.value,1)
    w.value{kk}= w.value{kk} + u.value{ll,kk};
        end
    end
elseif k==2
    w.sm = [size(u.value,1),1];
    w.value = cell(size(u.value,1),1);
    for kk=1:size(u.value,1)
    w.value{kk} = u.value{kk,1};    
        for ll=1:size(u.value,2)
    w.value{kk}= w.value{kk} + u.value{kk,ll};
        end
    end
end

end
end