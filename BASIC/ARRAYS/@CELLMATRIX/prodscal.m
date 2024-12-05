function an = prodscal(a,b,M,M2)


    if isa(a,'double')
     if nargin==3 && ~isempty(M)
         if ~isa(M,'double')
             error('pas programme')
         end        
         a=M'*a;
     end
     if nargin==4 && ~isempty(M2)
         if ~isa(M2,'double')
             error('pas programme')
         end        
         a=a*M2;
     end
     an=b;
     an.value = a(:)'*b.value ;
     an.s=[1,1];
     an.sm=b.sm;
    elseif isa(b,'double')
     if nargin==3 && ~isempty(M)
         if isa(M,'MULTIMATRIX')
             error('pas programme')
         end        
         b=M*b;
     end
     if nargin==4 && ~isempty(M2)
         if ~isa(M2,'double')
             error('pas programme')
         end        
         b=b*M2;
     end
     
     an=a;
     an.value = b(:)'*a.value ;
     an.s=[1,1];
     an.sm=a.sm;  
    else
    if nargin==3 && ~isempty(M)
         if isa(M,'MULTIMATRIX')
             error('pas programme')
         end        
         b=M*b;
     end
     if nargin==4 && ~isempty(M2)
         if ~isa(M2,'double')
             error('pas programme')
         end        
         b=b*M2;
     end
    
    an = sum(sum(times(a,b))); 
    %an = a ; 
    %an.value = a.value'*b.value ;
    %an.value = an.value(:)';
    %an.s=[1,1];
    %an.sm=[a.sm,b.sm];
    end
    
    