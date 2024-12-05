function Z = expectmtimes(a,b,c)
% function Z = expectmtimes(a,b)
% calcul de E(a*b)
%
% function Z = expectmtimes(a,b,c)
% calcul de E(a*b*c)

switch nargin
    case 2
  if ~israndom(a)
      Z = a*expect(b);
  elseif ~israndom(b)
      Z = expect(a)*b;
  elseif isa(a,'PCMATRIX') && isa(b,'PCMATRIX')
      if all(size(a)==1)     
       Z = reshape(double(b)*double(a)',size(b));
      elseif all(size(b)==1)     
       Z = reshape(double(a)*double(b)',size(a));                  
      else
      Z = multisum(a.MULTIMATRIX*b.MULTIMATRIX);      
      end
  else
      error('pas defini')   
  end
    case 3
  if ~israndom(a)
      Z = a*expectmtimes(b,c);
  elseif ~israndom(b)
      if isa(b,'double')
      Z = expectmtimes(a,mtimes(b,c));
      elseif isa(b,'MULTIMATRIX')
      Z = expectmtimes(a,mtimes(b{1},c));
      sz = size(Z);
      Z = Z(:);
      for i=1:prod(sizem(b))    
      Zi = expectmtimes(a,mtimes(b{i},c));    
      Z = [Z,Zi(:)];
      end
      Z = MULTIMATRIX(Z,sz,sizem(b));
      end
  elseif ~israndom(c)
      Z = expectmtimes(a,b)*c;     
  elseif isa(a,'PCMATRIX') && isa(b,'PCMATRIX') && isa(c,'PCMATRIX')
      if all(size(a)==1)     
       Z = expectmtimes(a*b,c);
      elseif all(size(b)==1) || all(size(c)==1)
      
       Z = expectmtimes(a,b*c);
     
      else
          if prod(size(a))<=prod(size(c))
             Z = expectmtimes(a*b,c) ;   
          else
              Z = expectmtimes(a,b*c);    
          end
      end
  else
      error('pas defini')   
  end        
    
        
end

