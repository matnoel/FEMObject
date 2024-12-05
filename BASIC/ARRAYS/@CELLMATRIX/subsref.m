function u = subsref(u,s)
% function u = subsref(u,s)

if length(s)==1 && (strcmp(s.type,'{}'))
    
    if length(s.subs)==2
        
        if isa(s.subs{1},'char') && strcmp(s.subs{1},':')
            s.subs{1} = 1:u.sm(1);
        end
        if isa(s.subs{2},'char') && strcmp(s.subs{2},':')
            s.subs{2} = 1:u.sm(2);
        end
        
        if length(s.subs{1})==1 && length(s.subs{2})==1
            I = sub2ind(u.sm,s.subs{1},s.subs{2});
            u = reshape(u.value(:,I),u.s);
        else
            u.sm = [length(s.subs{1}),length(s.subs{2})];
            s.subs{1} = repmat(s.subs{1}(:),1,length(s.subs{2}));
            s.subs{2} = repmat(s.subs{2}(:)',length(s.subs{1}),1);
            I = sub2ind(u.sm,s.subs{1},s.subs{2});
            u.value = u.value(:,I);
            
        end
        
        
    elseif length(s.subs)==1
        
        if length(s.subs{1})==1 && ~isa(s.subs{1},'char')
            u = reshape(u.value(:,s.subs{1}),u.s);
        else
            u.value = u.value(:,s.subs{1});
            if u.sm(1)==1
                u.sm = [1,size(u.value,2)];
            else
                u.sm = [size(u.value,2),1];
            end
        end
        
    end
    
elseif length(s)==2 && (strcmp(s(1).type,'{}')) && (strcmp(s(2).type,'()'))
    
    if length(s(1).subs)==1
        u.value = u.value(:,s(1).subs{1});
        u.sm = [size(u.value,2),1];
        if size(u.value,2)==1
            u = reshape(u.value,u.s);
            u = subsref(u,s(2));
        else
            u = subsref(u,s(2));
        end
    else
        error('pas programme')
    end
    
elseif length(s)==1 && strcmp(s.type,'()')
    
    if length(s.subs)==1 && length(size(u))==2
        
        if isa(s.subs{1},'char') && strcmp(s.subs{1},':')
            s.subs{1} = 1:numel(u);
            u.s = [numel(u),1];
        elseif u.s(1)==1
            u.s = [1,length(s.subs{1})];
        else
            u.s = [length(s.subs{1}),1];
        end
        
        u.value = u.value(s.subs{1},:);
        
    elseif length(s.subs)==2 && length(size(u))==2
        for i=1:2
            if isa(s.subs{i},'char') && strcmp(s.subs{i},':')
                s.subs{i} = 1:u.s(i);
            end
        end
        n = [length(s.subs{1}),length(s.subs{2})];
        I = repmat(s.subs{1}(:),n(2),1);
        J = repmat(s.subs{2}(:)',n(1),1);J=J(:);
        S = sub2ind(u.s,I,J);
        u.value = u.value(S,:);
        u.s = n ;
        
    elseif length(s.subs)==3 && length(size(u))==2
        error('pas programme')
    elseif length(s.subs)>=2 && length(size(u))==3
        if length(s.subs)==2
            s.subs{3} = 1:u.s(3);
        end
        for i=1:3
            if isa(s.subs{i},'char') && strcmp(s.subs{i},':')
                s.subs{i} = 1:u.s(i);
            end
        end
        n = [length(s.subs{1}),length(s.subs{2}),length(s.subs{3})];
        I = repmat(s.subs{1}(:),n(2)*n(3),1);
        J = repmat(s.subs{2}(:)',n(1),n(3));J=J(:);
        K = repmat(s.subs{3}(:)',n(1)*n(2),1);K=K(:);
        S = sub2ind(u.s,I,J,K);
        u.value = u.value(S,:);
        u.s = n ;
        
        
    else
        error('bad subsref')
        
    end
    
    
else
    error('bad subsref')
end


