function u=max(u,a,k)

if isa(u.value,'cell')
if nargin==1
 if ~all(u.s==1)
    for i=1:numel(u.value)
       u.value{i}=max(u.value{i}); 
    end     
    u.s = size(u.value{1});
 else
    u = [u.value{:}];
    u = max(u);
 end
     
elseif (nargin==3 & isempty(a))

    if k<=length(s)
    for i=1:numel(u.value)
       u.value{i}=max(u.value{i},[],k);
    end     
    u.s = size(u.value{1});
    else
      error('pas programme')  
    end
elseif isa(u,'MULTIMATRIX') & isa(a,'double')
    
    for i=1:numel(u.value)
       u.value{i}=max(u.value{i},a); 
    end     
elseif isa(u,'MULTIMATRIX') & isa(a,'MULTIMATRIX')
    if all(u.sm==a.sm)
     for i=1:numel(u.value)
       u.value{i}=max(u.value{i},a.value{i}); 
    end        
    else
        error('pas les meme multidimensions')
    end
    
end
    
else
if (nargin==3 & isempty(a))| (nargin==1) 
    
if nargin==1 
    a=[];
    r=find([u.s,u.sm]>1);
    if isempty(r)
     k=[];
    else
     k = r(1);
    end
end

s=u.s;
sm=u.sm;
%switch k
%    case 0
%    u.value = max(u.value,[],1);
%    u.s = [1,1]; 
%    case 1
%     u.value = max(reshape(u.value,s(1),s(2)*prod(sm),1)); 
%     u.s = [1,s(2)];
%     u.value = reshape(u.value,s(2),prod(sm));
%    case 2
%     u.value = max(reshape(u.value.',s(1)*prod(sm),s(2),2)); 
%     u.s = [s(1),1];
%     u.value = reshape(u.value,s(1),prod(sm)).';     
%    otherwise
     u.value=reshape(full(u.value),[s sm]);
     
     u.value=max(u.value,[],k);
     
     if k<=length(s)
     s(k)=1;
     else
     sm(k-length(s))=1;
     end    
     
    
     if all(u.s)==1 && prod(sm)==1
     u = u.value ;  
     else
     u.value = reshape(u.value,prod(s),prod(sm));
     u.s=s;
     u.sm=sm;
     end
%end


elseif isa(u,'MULTIMATRIX') & isa(a,'double')
if all(size(a)==1)    
   u.value = max(u.value,a);
elseif all(size(a)==size(u))   
   u.value = max(u.value,repmat(a(:),1,prod(u.sm)));
else
error('MULTIMATRIX dimensions must agree')
end
elseif isa(u,'MULTIMATRIX') & isa(a,'MULTIMATRIX')
if all(size(u)==size(a))
   u.value = max(u.value,a.value);
else
error('MULTIMATRIX dimensions must agree')
end
        
else
    
    error('bad arguments')
    
end
end