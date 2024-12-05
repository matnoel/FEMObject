function u=getmultimatrix(u,k,l)

switch nargin
    case 2
if ~isa(k,'char') && ~strcmp(k,':') 
if all(size(k)>1)
   error('rentrer un vecteur') 
end    
    
 if isa(u.value,'cell')
     u.value=u.value(k);     
 else
     u.value=u.value(:,k);     
 end
 if u.sm(1)==1
 u.sm = [1,length(k)];
 else
 u.sm = [length(k),1];
 end 
end

    case 3
 if ~isa(k,'char') && ~strcmp(k,':')  && ~isa(l,'char') && ~strcmp(l,':')    
  if all(size(k)>1) ||  all(size(l)>1)
   error('rentrer un vecteur') 
  end   
 K = repmat(k(:),1,length(l));
 L = repmat(l(:)',length(k),1);
 I=sub2ind(u.sm,K(:),L(:));
 if isa(u.value,'cell')
     u.value=u.value(I);     
 else
     u.value=u.value(:,I);     
 end
 u.sm = [length(k),length(l)];
  
 end
end