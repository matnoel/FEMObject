function Z = expecttimes(a,b,c)


switch nargin
    case 2
  if ~israndom(a) & ~israndom(b)
      Z = times(a,b);    
  elseif ~israndom(a)
      Z = times(a,expect(b));
  elseif ~israndom(b)
      Z = times(expect(a),b);
  else 
      error('pas possible')
  end
    case 3
  if ~israndom(a)
      Z = times(a,expecttimes(b,c));
  elseif ~israndom(c)
      Z = expecttimes(a,b)*c;     
  elseif ~israndom(b)
      Z = expecttimes(a,times(b,c));
  else
  error('pas possible')
  end        
    
        
end

