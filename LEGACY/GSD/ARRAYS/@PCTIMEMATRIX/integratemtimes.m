function Z = integratemtimes(a,b,c)


switch nargin
    case 2
if ~istime(b)
 Z = mtimes(integrate(a),b);   
elseif ~istime(a)
 Z = mtimes(a,integrate(b));       
elseif isa(a,'TIMEMATRIX')
  if isa(b.value,'PCRADIALMATRIX')
   B = TIMEMATRIX(getV(b.value),b.TIMEMODEL,b.s);

   Z = integratemtimes(a,B);
   Z = PCRADIALMATRIX(Z,size(Z),getL(b.value));
   Z = setmasse(Z,getmasse(b.value));
   
   %B = TIMEMATRIX(getV(mat2cell(b.value)),b.TIMEMODEL,b.s);
   %Z = integratemtimes(a,B);
   %Z = PCRADIALMATRIX(Z,size(Z),getL(b.value));
   %Z = setmasse(Z,getmasse(b.value));
  elseif isa(b.value,'PCMATRIX')
   B = TIMEMATRIX(getmultimatrix(b.value),b.TIMEMODEL,b.s);
   Z = integratemtimes(a,B);
   Z = PCMATRIX(Z,size(Z),getPC(b));   
   
  end
  
elseif isa(b,'TIMEMATRIX')
  if isa(a.value,'PCRADIALMATRIX')
   A = TIMEMATRIX(getV(a.value),a.TIMEMODEL,a.s);
   Z = integratemtimes(A,b);
   Z = PCRADIALMATRIX(Z,size(Z),getL(a.value));
   Z = setmasse(Z,getmasse(a.value));
  elseif isa(b.value,'PCMATRIX')
   A = TIMEMATRIX(getmultimatrix(a.value),a.TIMEMODEL,a.s);
   Z = integratemtimes(A,b);
   Z = PCMATRIX(Z,size(Z),getPC(a));   
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
    error('pas programme')
end

    
    
end





