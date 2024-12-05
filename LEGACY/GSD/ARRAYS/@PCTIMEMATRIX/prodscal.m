function an=prodscal(a,b,A)

switch nargin
    case 3
if ~istime(A) & ~israndom(A)
    an = prodscal(a,mtimes(A,b));
elseif ~istime(A)
    if ~israndom(a)
    an = prodscal(b,mtimes(A,a));    
    else
    an = prodscal(a,mtimes(A,b));        
    end
else
   error('pas programme') 
end
    
    case 2
if ~israndom(a)
    an = prodscal(a,expect(b));
elseif ~israndom(b)
    an = prodscal(expect(a),b);
elseif isa(a,'PCTIMEMATRIX') & isa(b,'PCTIMEMATRIX')
    if isa(a.value,'cell') & ~isa(b.value,'cell')
    b = mat2cell(b);
    end
    if isa(b.value,'cell') & ~isa(a.value,'cell')
    a = mat2cell(a);
    end

    an = prodscal(a.value,b.value,[],getMmatrix(a));
    
elseif ~istime(a)  
    an = prodscal(a,integrate(b));
elseif ~istime(b)
    an = prodscal(integrate(a),b);
else
    error('pas programme')
end
 
end
   
an=full(an);
    

%   an = b.MULTIMATRIX(:)'*a.MULTIMATRIX(:);
%   s = size(an);
%   an = reshape(full(sum(double(an),2)),s); 

