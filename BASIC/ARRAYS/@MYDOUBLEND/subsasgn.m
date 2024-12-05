function u = subsasgn(u,s,b)
% function u = subsasgn(u,s,b)

switch s.type
    case {'()'}
        v = double(u);
        b = double(b);
        if length(s.subs)==1
            v = reshape3D(v);
            v = reshape(v,[size(v,1)*size(v,2),1,size(v,3)]);
            b = samesize(b,double(subsref(u,s)));
            b = reshape3D(b);
            b = reshape(b,[size(b,1)*size(b,2),1,size(b,3)]);
            v(s.subs{1},1,:) = b;
            v = reshape(v,size(u));
        elseif length(s.subs)==2
            v = reshape3D(v);
            b = samesize(b,double(subsref(u,s)));
            b = reshape3D(b);
            v(s.subs{1},s.subs{2},:) = b;
            v = reshape(v,size(u));
        else
            b = samesize(b,double(subsref(u,s)));
            v = subsasgn(v,s,b);
        end
        u.double = v;
end
