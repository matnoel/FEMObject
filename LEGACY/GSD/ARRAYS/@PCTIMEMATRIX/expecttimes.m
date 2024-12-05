function w = expecttimes(a,b,c)

switch nargin
    case 2
  if ~israndom(a)
      Z = a.*expect(b);
  elseif ~israndom(b)
      Z = expect(a).*b;
  else
    if isa(a,'TIMEMATRIX') & isa(b,'TIMEMATRIX')
        
    elseif isa(b,'TIMEMATRIX')
    if all(size(a)==1)
      Z = b ; 
      Z.value = expecttimes(a,b.value);
    else
        error('pas programme')
    end
    elseif isa(a,'TIMEMATRIX')        
    if all(size(b)==1)
      Z = a ; 
      Z.value = expecttimes(a.value,b);
    else
        error('pas programme')
    end       
    end

  end
    case 3
  if ~israndom(a)
      Z = a*expecttimes(b,c);
  elseif ~israndom(b)
      Z = b*expecttimes(a,c);
  elseif ~israndom(c)
      Z = expecttimes(a,b).*c;     
  elseif israndom(a) & israndom(b) & israndom(c)
      if all(size(a)==1)     
       Z = expecttimes(a.*b,c);
      elseif all(size(b)==1) | all(size(c)==1)
       Z = expecttimes(a,b.*c);     
      else
          if prod(size(a))<=prod(size(c))
             Z = expecttimes(a.*b,c) ;   
          else
              Z = expecttimes(a,b.*c);    
          end
      end
  else
      error('pas defini')   
  end        
    
        
end

