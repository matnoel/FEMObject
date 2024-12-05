function an=prodscal(a,b,A)

if nargin==3   
if ~istime(A) && ~israndom(A)
    an = prodscal(a,A*b); 
else
    error('pas programme')
end

elseif nargin==2
if isa(a,'TIMEMATRIX') & ~isa(b,'TIMEMATRIX')    
    an = prodscal(integrate(a),b);
elseif ~isa(a,'TIMEMATRIX') & isa(b,'TIMEMATRIX')      
    an = prodscal(a,integrate(b));
elseif isa(a,'TIMEMATRIX') & isa(b,'TIMEMATRIX')
     %Mt = getMmatrix(a);   
     %a.value = a.value*Mt;
     
     an = prodscal(a.value,b.value,[],getMmatrix(a));
else
    error('pas programme')
end
  
end
   

 an=full(an);
  
%   an = b.MULTIMATRIX(:)'*a.MULTIMATRIX(:);
%   s = size(an);
%   an = reshape(full(sum(double(an),2)),s); 

