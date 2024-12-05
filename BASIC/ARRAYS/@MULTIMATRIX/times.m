function w=times(u,v,varargin)

if nargin>=3
    masse=varargin{1};
end

if nargin>=4
    n=varargin{2};
end



if isa(u,'MULTIMATRIX') & isa(v,'MULTIMATRIX') & all(u.sm==v.sm)
    
    if isa(u.value,'double') & isa(v.value,'cell')
    w = times(mat2cell(u),v,varargin{:});
        
    elseif isa(u.value,'cell') & isa(v.value,'double')
    w = times(u,mat2cell(v),varargin{:});
        
    elseif isa(u.value,'cell') & isa(v.value,'cell')
        if nargin==2        
           w=u; 
           
         for k=1:numel(u.value)
        w.value{k} = sparsetimes(u.value{k},v.value{k});
        w.s = size(w.value{k});
         end
      
         
         
        else
            error('pas programme')
            
        end
        
    else
    s = max(u.s,v.s);
    
    if nargin==2        
    p=size(u.value,2);
    w=u;
    w.s=s;
    for k=1:p
    w.value(:,k) = u.value(:,k).*v.value(:,k);
    end
    
    else
    value = [] ;
    i=[];
    j=[];
    if isa(masse,'MULTIMATRIX') & (nargin==3 | n==1)
    p=size(v.value,2);
    multivtemp = v.value*masse ;
    for k=1:size(u.value,2)
    vtemp = getmatrix(multivtemp,k);
    
    uvtemp = repmat(u.value(:,k),1,p).*vtemp;
    [ik,jk]=find(uvtemp);
    i=[i;ik];
    j=[j;jk];
    Ik = find(uvtemp);
    uvtemp = uvtemp(Ik);
    value=[value;uvtemp(:)];
    end
    %value = sparse(i,j,value,prod(s),p);
    %w = u ; 
    %w.value = value ; 
    %w.s = s;
    
    elseif isa(masse,'MULTIMATRIX') & n==2
    p=size(u.value,2);
    multivtemp = u.value*masse ;
    
    for k=1:size(v.value,2)    
    utemp = getmatrix(multivtemp,k);
    uvtemp = utemp.*repmat(v.value(:,k),1,p);
    [ik,jk]=find(uvtemp);
    
    i=[i;ik];
    j=[j;jk];
    
    Ik = find(uvtemp);
    uvtemp = uvtemp(Ik);
    value=[value;uvtemp(:)];

    end    
    end
    value = sparse(i,j,value,prod(s),p);
    w = u ; 
    w.value = value ; 
    w.s = s;
    end
    
    end
    
elseif isa(u,'MULTIMATRIX') & isa(v,'double')
    
    if isa(u.value,'cell')
        w=u;
        for k=1:numel(u.value)
        w.value{k} = sparsetimes(u.value{k},v);
        w.s = size(w.value{k});
        end
        
        
    else
    if all(size(v)==1)
    s=u.s;
    value = u.value.*v;
    else    
    p=size(u.value,2);
    s = u.s;
    value = u.value.*repmat(v(:),1,p);
    end
    
    w = u ; 
    w.value = value ; 
    w.s = s;
    end
    
elseif isa(v,'MULTIMATRIX') & isa(u,'double')
    if isa(v.value,'cell')
        w=v;
        
        for k=1:numel(v.value)
        w.value{k} = sparsetimes(u,v.value{k});
        w.s = size(w.value{k});
        end
      
        
    else
    if all(size(u)==1)
    s=v.s;
    value = u*v.value;
    else    
    p=size(v.value,2);
    s = v.s;
    value = repmat(u(:),1,p).*v.value;
    end

    w = v ; 
    w.value = value ; 
    w.s = s;
    end
elseif isa(u,'MULTIMATRIX') & isa(v,'MULTIMATRIX') & ~all(u.sm==v.sm)
if length(u.sm)>2 | length(v.sm)>2
    error('pas programme')
end
if nargin>2
   error('pas programme') 
end

    if isa(u.value,'double') & isa(v.value,'cell')
    w = times(mat2cell(u),v,varargin{:});
        
    elseif isa(u.value,'cell') & isa(v.value,'double')
    w = times(u,mat2cell(v),varargin{:});
        
    elseif isa(u.value,'cell') & isa(v.value,'cell')
if all(u.sm==1)
    u = reshapem(u,v.sm);
elseif all(v.sm==1)
    v = reshapem(v,u.sm); 
elseif u.sm(1)==1 & v.sm(2)==1
    u.value = repmat(u.value,v.sm(1),1);
    u.sm = size(u.value);
    v.value = repmat(v.value,1,u.sm(2));
    v.sm = size(v.value);
elseif u.sm(2)==1 & v.sm(1)==1
    u.value = repmat(u.value,1,v.sm(2));
    u.sm = size(u.value);
    v.value = repmat(v.value,u.sm(1),1);
    v.sm = size(v.value);
else
    error('pas programme')
end

    w = times(u,v);

    else

    
if all(u.sm==1)
    w = reshape(u.value,u.s)*v
elseif all(v.sm==1)
w = u*reshape(v.value,v.s);
elseif u.sm(1)==1 & v.sm(2)==1
v.value = repmat(v.value,1,u.sm(2));
u.value = repmat(u.value,1,v.sm(1));
sm = [v.sm(1),u.sm(2)];
u.sm=sm;
v.sm=sm;
u = multitranspose(u);
w = times(u,v);
    
elseif u.sm(2)==1 & v.sm(1)==1

v.value = repmat(v.value,1,u.sm(1));
u.value = repmat(u.value,1,v.sm(2));
sm = [u.sm(1),v.sm(2)];
u.sm=sm;
v.sm=sm;
v=multitranspose(v);
w = times(u,v);
 
else
    error('pas defini')
end

    end
    
else
       error('times not defined') 
end


