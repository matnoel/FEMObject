function an = prodscal(a,b,varargin)
if nargin>=3
    M = varargin{1};
end
if nargin>=4
    M2 = varargin{2};
end



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
     if isa(b.value,'cell')
       for k=1:numel(b.value)
        an.value{k} = a(:)'*b.value{k} ;   
       end
     else
     an.value = a(:)'*b.value ;
     end
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
     if isa(a.value,'cell')
       for k=1:numel(a.value)
        an.value{k} = b(:)'*a.value{k} ;   
       end
     else
     an.value = b(:)'*a.value ;
     end
     an.s=[1,1];
     an.sm=a.sm;  
 
elseif isa(a,'MULTIMATRIX') & isa(b,'MULTIMATRIX')

    if isa(a.value,'cell') & isa(b.value,'double')
        an = prodscal(a,mat2cell(b),varargin{:});
    elseif isa(b.value,'cell') & ~isa(a.value,'cell')
       an = prodscal(mat2cell(a),b,varargin{:});
    elseif isa(b.value,'cell') & isa(a.value,'cell')
        if ~all(a.sm==b.sm)
            error('pas la meme multidim')
        end
        an = a;
        
        for k=1:numel(a.value)
        an.value{k} = prodscal(a.value{k},b.value{k},varargin{:}) ;   
        end
       an.s = [1,1];
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
    end
    
    %an = a ;
    %an.value = a.value'*b.value ;
    %an.value = an.value(:)';
    %an.s=[1,1];
    %an.sm=[a.sm,b.sm];
    
else
    error('erreur pas programme')
    
end
    
    