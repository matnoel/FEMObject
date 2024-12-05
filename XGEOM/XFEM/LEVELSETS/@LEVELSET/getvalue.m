function value = getvalue(ls,listenode)
% function value = getvalue(ls,listenode)

if ~iseval(ls)
    error('value definie pour une LEVELSET')
else
    if nargin==1
        value=ls.value;
        if isa(value,'MULTIMATRIX') && length(value)==1 && ~israndom(value)
            value = double(value) ;
        end
    else
        if israndom(ls.value)
            error('a voir')
        end
        ls.value = double(ls.value);
        if any(size(listenode)==1)
            listenode=listenode(:);
        end
        pos=listenode;
        lsx = ls.value;
        value = zeros([size(pos),size(lsx,2)]);
        coco = repmat(pos(:),1,size(lsx,2)) + int32(repmat(size(lsx,1),numel(pos),1)*[0:size(lsx,2)-1]);
        value(:) = lsx(coco(:));
        value = permute(value,[1,4,2,3]);
        value = MYDOUBLEND(value);
        
    end
end