function Z = integratetimes(a,b,c)


switch nargin
    case 2
if ~istime(b)
 Z = times(integrate(a),b);   
elseif ~istime(a)
 Z = times(a,integrate(b));       
elseif isa(a,'TIMEMATRIX') && isa(b,'TIMEMATRIX')
  if isa(a.value,'double') && isa(b.value,'double')
   M=getMmatrix(a);
   
   if all(size(a)==1) && all(size(b)==1) 
   Z = a.value(:)'*M*b.value(:);
   elseif all(size(a)==1)
   Z = b.value*M*a.value(:);    
   elseif all(size(b)==1)
   Z = a.value*M*b.value(:);     
   else
   A = a.value*M;
   A = MULTIMATRIX(A,a.s,[length(a.TIMEMODEL),1]);
   B = MULTIMATRIX(b.value,b.s,[length(b.TIMEMODEL),1]);
   Z = multisum(A.*B);
    
   end
   
  else    
    if isa(a.value,'MULTIMATRIX') && isa(b.value,'MULTIMATRIX') && ~all(sizem(a.value)==sizem(b.value))
error('les MULTIMATRIX doivent avoir les memes multi-dimensions')
    end
  A = MULTIMATRIX(double(a.value),a.s,[length(a.TIMEMODEL),prod(sizem(a.value))]);
  B = MULTIMATRIX(double(b.value),b.s,[length(b.TIMEMODEL),prod(sizem(b.value))]);
 
  Z = multisum(A.*B,1);
  if isa(a.value,'MULTIMATRIX') 
   Z = reshapem(Z,sizem(a.value));
  elseif isa(b.value,'MULTIMATRIX') 
   Z = reshapem(Z,sizem(b.value));
  else
   Z = double(Z);   
  end
  end
else
    error('pas programme')
end    
    case 3
if ~istime(a) && ~istime(b)
    Z = times(times(a,b),integrate(c));
elseif ~istime(b) && ~istime(c)
    Z = times(integrate(a),times(b,c));
elseif ~istime(a) && ~istime(c)
    Z = times(a,times(integrate(b),c));
elseif ~istime(a)
  Z = times(a,integratetimes(b,c));  
elseif ~istime(c)
  Z = times(integratetimes(a,b),c); 
elseif ~istime(b)
  Z = integratetimes(a,times(b,c));    
else
  Z = integratetimes(a,times(b,c));  
end
end

