function Z = expectmtimes(a,b,c)

switch nargin
    case 2
  if ~israndom(a)
      Z = a*expect(b);
  elseif ~israndom(b)
      Z = expect(a)*b;
  else
    if isa(a,'PCTIMEMATRIX') && isa(b,'PCTIMEMATRIX')
        error('pas programme')
    elseif isa(b,'PCTIMEMATRIX')
    if all(size(a)==1)
      Z = expectmtimes(a,b.value);
      Z = TIMEMATRIX(Z,b.TIMEMODEL,b.s);  
    else
        error('pas programme')
    end
    elseif isa(a,'PCTIMEMATRIX')        
    if all(size(b)==1)
      Z = expectmtimes(a.value,b);
      Z = TIMEMATRIX(Z,a.TIMEMODEL,a.s);
    else
        error('pas programme')
    end       
    end

  end
    case 3
  if ~israndom(a)
      Z = a*expectmtimes(b,c);
  elseif ~israndom(c)
      Z = expectmtimes(a,b)*c;     
  elseif ~israndom(b)
      if isa(b,'double')
      Z = expectmtimes(a,mtimes(b,c));
      elseif isa(b,'TIMEMATRIX') & ~isa(c,'PCTIMEMATRIX')
      Z = expectmtimes(a,mtimes(b,c));    
      elseif isa(b,'TIMEMATRIX') & ~isa(a,'PCTIMEMATRIX')
      Z = expectmtimes(mtimes(a,b),c); 
      else
          error('pas programme')
      end
  elseif israndom(a) & israndom(b) & israndom(c)
      if all(size(a)==1) || prod(size(a))<=prod(size(c))
       Z = expectmtimes(a*b,c);
      else
       Z = expectmtimes(a,b*c);     

      end
  else
      error('pas defini')   
  end        
    
        
end

