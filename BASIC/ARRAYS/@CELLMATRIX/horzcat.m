function w = horzcat(u,v,varargin)
% function w = horzcat(u,v,varargin)

if nargin==1
    w = u;
else
    if isa(u,'double')
        p = size(v.value,2);
        s = size(u);
        u = repmat(u(:),1,p);
        u = MULTIMATRIX(u,s,sizem(v));
        w = horzcat(u,v);
        
    elseif isa(v,'double')
        p = size(u.value,2);
        s = size(v);
        v = repmat(v(:),1,p);
        v = MULTIMATRIX(v,s,sizem(u));
        w = horzcat(u,v);
        
    elseif isa(u,'MULTIMATRIX') && isa(v,'MULTIMATRIX')
        w = u;
        if u.s(1)==v.s(1)
            w.s = [u.s(1),u.s(2)+v.s(2)];
            w.value = [u.value;v.value];
        else
            error('pas la bonne dimension')
        end
        
%     else
%         l = [u.s(1),v.s(1)];
%         c = [u.s(2),v.s(2)];
%         p = [size(u.value,2),size(v.value,2)];
%         value = {u.',v.'};
%         for k=1:2
%             value{k} = reshape(value{k}.value,c(k),l(k)*p(k));
%         end
%         
%         
%         value = vertcat(value{:});
%         value = reshape(value,[sum(c)*l(1)],p(1));
%         w.s = [sum(c) , l(1)] ;
%         w.value = value ;
%         w = w.';
%     end
        
    else
        error('pas defini')
    end
    
    if nargin>2
        w = horzcat(w,varargin{:});
    end
    
end
