function w=vertcat(u,v,varargin)

if nargin==1
   w=u;
else
    
    if isa(u,'double')
        if isa(v.value,'cell')
     w=v;
     for k=1:numel(w.value)
     w.value{k}=vertcat(u,w.value{k});    
     end
     w.s = size(w.value{1});
        else
            p = size(v.value,2);
     s=size(u);
     u = repmat(u(:),1,p);
     u = MULTIMATRIX(u,s,sizem(v));
     w = vertcat(u,v);
        end
    elseif isa(v,'double')
     if isa(u.value,'cell')
     w=u;
     for k=1:numel(w.value)
     w.value{k}=vertcat(w.value{k},v);    
     end
     w.s = size(w.value{1});
     else
         p = size(u.value,2);
     s=size(v);
     v = repmat(v(:),1,p);
     v = MULTIMATRIX(v,s,sizem(u));
     w = vertcat(u,v);
     end
    elseif isa(u,'MULTIMATRIX') & isa(v,'MULTIMATRIX')
      if isa(u.value,'cell') &  ~isa(v.value,'cell')
          v = mat2cell(v);
      elseif ~isa(u.value,'cell') &  isa(v.value,'cell')
          u = mat2cell(u);
      end
      
      if isa(u.value,'cell') &  isa(v.value,'cell')
        if all(u.sm==v.sm)
            w=u;
            for k=1:numel(w.value)
            w.value{k}=vertcat(u.value{k},v.value{k});
            end
            w.s = size(w.value{1});
         else
             error('pas les memes multidimensions')
        end
         
        elseif isa(u.value,'double') &  isa(v.value,'double')   
          w=u;
     if u.s(2)==1 && v.s(2)==1
     w.s=[u.s(1)+v.s(1),u.s(2)];
     w.value = [u.value;v.value];
     else
     l = [u.s(1),v.s(1)];
     c = [u.s(2),v.s(2)];
     p = [size(u.value,2),size(v.value,2)];
     value = {u.value,v.value};
     for k=1:2
     value{k}=reshape(value{k},l(k),c(k)*p(k));  
     end
     
     
    value = vertcat(value{:});
    value = reshape(value,[sum(l)*c(1)],p(1));
    w.s = [sum(l) , c(1)] ; 
    w.value = value ;
     end
     
     else
         error('horzcat pour 2 cell ou 2 double')
     end
    else
        error('pas defini')
    end
    
if nargin>2
  w = vertcat(w,varargin{:});
end
  
end



