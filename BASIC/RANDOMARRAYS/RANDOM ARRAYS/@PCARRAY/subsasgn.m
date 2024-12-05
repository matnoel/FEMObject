function u = subsasgn(u,s,v)
% function u = subsasgn(u,s,v)

if length(s)==1 && strcmp(s.type,'()') && length(s.subs)==ndims(u)
    
    u.MYDOUBLE = subsasgn(u.MYDOUBLE,s,v);
    
elseif length(s)==1 && strcmp(s.type,'()')
    
    nopc = setdiff(1:ndims(u),u.pcdim);
    rep = cell(1,ndims(u));
    rep(nopc) = s.subs;
    rep{u.pcdim} = ':';
    u.MYDOUBLE(rep{:}) = v ;
    
elseif length(s)==1 && strcmp(s.type,'{}') && length(s.subs)==1
    rep = cell(1,ndims(u));
    rep(:) = {':'};
    rep{u.pcdim} = s.subs{1};
    u.MYDOUBLE(rep{:}) = v;
    
else
    error('bad subsref')
end
