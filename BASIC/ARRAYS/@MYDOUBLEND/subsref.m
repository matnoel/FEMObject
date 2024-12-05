function u = subsref(u,s)
% function u = subsref(u,s)

switch s.type
    case {'()'}
        if length(s.subs)==1
            v = reshape3D(u);
            v = reshape(v,[size(v,1)*size(v,2),1,size(v,3)]);
            v = v(s.subs{1},1,:);
            v = reshape(v,[size2D(v),sizeND(u)]);
        elseif length(s.subs)==2
            v = reshape3D(u);
            v = v(s.subs{1},s.subs{2},:);
            v = reshape(v,[size2D(v),sizeND(u)]);
        else
            v = subsref(double(u),s);
        end
        u.double = v;
end

