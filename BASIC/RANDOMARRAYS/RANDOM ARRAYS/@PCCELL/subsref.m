function v = subsref(u,s)
% function v = subsref(u,s)

if length(s)==1 && (strcmp(s.type,'{}'))
    if length(s.subs)==1
        v = u.value{s.subs{1}};
    end
elseif length(s)==1 && strcmp(s.type,'.')
    v = getfield(struct(u),s.subs) ;
elseif length(s)==1 && strcmp(s.type,'()')
    si = size(u);
    v = [u.value{:}];
    if length(s.subs)==2
        ind = sub2ind(si,s.subs{1},s.subs{2});
    elseif length(s.subs)==1
        ind = sub2ind(si,s.subs{1});
    end
    v = reshape(v,prod(si),getP(u)+1);
    v = PCARRAY(v(ind,:),getPC(u));
else
    error('bad subsref')
end