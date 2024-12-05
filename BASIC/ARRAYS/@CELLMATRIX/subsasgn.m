function u = subsasgn(u,s,v)
% function u = subsasgn(u,s,v)

if length(s)==1 && (strcmp(s.type,'{}'))
    if length(s.subs)==1
        u.value(:,s.subs{1}) = reshape(v,prod(size(v)),1);
    elseif length(s.subs)==2
        rep = sub2ind(sizem(u),s.subs{1},s.subs{2});
        u.value(:,rep) = reshape(v,prod(size(v)),1);
    else
        error('pas programme ')
    end
    
elseif length(s)==1 && strcmp(s.type,'()')
    
    if length(s.subs)==1
        if isa(s.subs{1},'char') && strcmp(s.subs{1},':')
            s.subs{1}=1:numel(u);
        end
        ind = s.subs{1}(:);
        if any(ind>prod(size(u)))
            error('index exceeds matrix dimension')
        end
        p=length(u);
        v=double(v);
        if numel(v)==numel(ind)
            v=repmat(v(:),1,p);
        end
        v=reshape(v,numel(ind),p);
        u.value(ind,:) = v ;
    elseif length(s.subs)>=2
        for k=1:2
            if isa(s.subs{k},'char') && strcmp(s.subs{k},':')
                s.subs{k} = 1:size(u,k);
            end
        end
        if any(s.subs{1}(:)>size(u,1)) || any(s.subs{2}(:)>size(u,2))
            error('index exceeds matrix dimension')
        end
        
        ind1=repmat(s.subs{1}(:),1,numel(s.subs{2}));
        ind2=repmat(s.subs{2}(:)',numel(s.subs{1}),1);
        ind = sub2ind(size(u),ind1(:),ind2(:));
        p=length(u);
        v=double(v);
        if numel(v)==numel(ind)
            v=repmat(v(:),1,p);
        end
        v=reshape(v,numel(ind),p);
        u.value(ind,:) = v ;
    end
    
else
    error('bad subsref')
end