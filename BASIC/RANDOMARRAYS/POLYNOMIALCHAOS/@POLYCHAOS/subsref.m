function v = subsref(u,s)
% function v = subsref(pc,s)
% pc : POLYCHAOS in dimension M
% x : array of size n*M
% pc(x) evaluate all basis functions (P+1) on n points defined by x
% -> v array of size n*(P+1)
% pc(i,x) evaluate the basis functions listed in i (in 0,...,P)
% on n points defined by x
% -> v array of size n * length(i)

if length(s)==1 && strcmp(s.type,'()')
    
    switch length(s.subs)
        case 2
            v = polyval(u,s.subs{2},s.subs{1});
        case 1
            v = polyval(u,s.subs{1});
    end
    
elseif length(s)==1 && strcmp(s.type,'.')
    v = getfield(struct(u),s.subs) ;
else
    error('bad subsref')
end

