function v = subsref(u,s)
% function v = subsref(u,s)

if length(s)==1 && strcmp(s.type,'()') && length(s.subs)==ndims(u)
    if strcmp(s.subs{u.pcdim},':')
        v = u;
        v.MYDOUBLE = subsref(u.MYDOUBLE,s);
    else
        v = subsref(double(u.MYDOUBLE),s);
    end
elseif length(s)==1 && strcmp(s.type,'()')
    
    nopc = setdiff(1:ndims(u),u.pcdim);
    rep = cell(1,ndims(u));
    rep(nopc) = s.subs;
    rep{u.pcdim} = ':';
    v = u;
    v.MYDOUBLE = u.MYDOUBLE(rep{:});
    
elseif length(s)==1 && strcmp(s.type,'{}') && length(s.subs)==1
    
    rep = cell(1,ndims(u));
    rep(:) = {':'};
    rep{u.pcdim} = s.subs{1};
    v = double(u.MYDOUBLE);
    v = v(rep{:});
    
else
    error('bad subsref')
end
