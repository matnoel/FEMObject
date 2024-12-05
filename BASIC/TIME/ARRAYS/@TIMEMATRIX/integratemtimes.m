function Z = integratemtimes(a,b,c)

switch nargin
    case 2
if ~istime(b)
 Z = mtimes(integrate(a),b);   
elseif ~istime(a)
 Z = mtimes(a,integrate(b));       
elseif isa(a,'TIMEMATRIX') && isa(b,'TIMEMATRIX')
 if isa(a.value,'double') && isa(b.value,'double')
   M=getMmatrix(a);
   if all(size(a)==1)
   A=a.value*M;
   Z = reshape(double(b.value)*A',size(b));
   elseif all(size(b)==1)
   B=b.value*M;
   Z = reshape(double(a.value)*B',size(a));
   else
   A = a.value*M;
   A = MULTIMATRIX(A,a.s,[length(a.TIMEMODEL),1]);
   B = MULTIMATRIX(b.value,b.s,[length(b.TIMEMODEL),1]);
   Z = multisum(A*B);
   end
   
 elseif isa(a.value,'cell') && isa(b.value,'cell')
     M = getMmatrix(a);
     A = getcell(a);
     B = getcell(b);
     Z = multicellmtimes(M,B(:));
     Z = cellmtimes(A(:)',Z(:));
     
         
 elseif isa(a.value,'cell') && ~isa(b.value,'cell')
     b = mat2cell(b);
     Z = integratemtimes(a,b);
 elseif ~isa(a.value,'cell') && isa(b.value,'cell')
     a = mat2cell(a);
     Z = integratemtimes(a,b);
 else
 
    
    if isa(a.value,'MULTIMATRIX') && isa(b.value,'MULTIMATRIX') && ~all(sizem(a.value)==sizem(b.value))
error('les MULTIMATRIX doivent avoir les memes multi-dimensions')
    end
  
  A = MULTIMATRIX(a.value,a.s,[length(a.TIMEMODEL),prod(sizem(a.value))]);
  B = MULTIMATRIX(b.value,b.s,[length(b.TIMEMODEL),prod(sizem(b.value))]);
  if iscell(A) && ~iscell(B)
     B=mat2cell(B);
  elseif ~iscell(A) && iscell(B)
     A=mat2cell(A);
  end
  
  Z = multisum(A*B,1);
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
    Z = mtimes(mtimes(a,b),integrate(c));
elseif ~istime(b) && ~istime(c)
    Z = mtimes(integrate(a),mtimes(b,c));
elseif ~istime(a) && ~istime(c)
    Z = mtimes(a,mtimes(integrate(b),c));
elseif ~istime(a)
  Z = mtimes(a,integratemtimes(b,c));  
elseif ~istime(c)
  Z = mtimes(integratemtimes(a,b),c); 
elseif ~istime(b)
  Z = integratemtimes(a,mtimes(b,c));    
else
  Z = integratemtimes(a,mtimes(b,c));  
end
end

