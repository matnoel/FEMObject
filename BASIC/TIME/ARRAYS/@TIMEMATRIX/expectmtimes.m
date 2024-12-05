function Z = expectmtimes(a,b,c)


switch nargin
    case 2
  if ~israndom(a) & ~israndom(b)
      Z = mtimes(a,b);    
  elseif ~israndom(a)
      Z = mtimes(a,expect(b));
  elseif ~israndom(b)
      Z = mtimes(expect(a),b);
  else
   error('pas possible')
  end
    case 3
  if ~israndom(a)
      Z = mtimes(a,expectmtimes(b,c));
  elseif ~israndom(c)
      Z = expectmtimes(a,b)*c;     
  elseif ~israndom(b)
      Z = expectmtimes(a,mtimes(b,c));
  else
   error('pas possible')
  end        
    
        
end

