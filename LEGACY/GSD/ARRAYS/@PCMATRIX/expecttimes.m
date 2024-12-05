function Z = expecttimes(a,b,c)
% function Z = expecttimes(a,b)
% calcul de E(a.*b)
%
% function Z = expecttimes(a,b,c)
% calcul de E(a.*b.*c)

switch nargin
    case 2
  if ~israndom(a)
      Z = a.*expect(b);
  elseif ~israndom(b)
      Z = expect(a).*b;
  elseif isa(a,'PCMATRIX') & isa(b,'PCMATRIX')
      if all(size(a)==1)     
       Z = reshape(double(b)*double(a)',size(b));
      elseif all(size(b)==1)
       Z = reshape(double(a)*double(b)',size(a));   
      else
      Z = multisum(a.MULTIMATRIX.*b.MULTIMATRIX);      
      end
  else
      error('pas defini')   
  end
    case 3
  if ~israndom(a)
      Z = a.*expecttimes(b,c);
  elseif ~israndom(b)
      Z = expecttimes(a,c).*b;
  elseif ~israndom(c)
      Z = expecttimes(a,b).*c;     
  elseif isa(a,'PCMATRIX') & isa(b,'PCMATRIX') & isa(c,'PCMATRIX')
      if all(size(a)==1)     
       Z = expecttimes(a.*b,c);
      elseif all(size(b)==1) | all(size(c)==1)
       Z = expecttimes(a,b.*c);     
      else
       if prod(size(a))<=prod(size(c)) & prod(size(b))<=prod(size(c))
             Z = expecttimes(a.*b,c) ;   
       elseif prod(size(c))<=prod(size(b)) & prod(size(a))<=prod(size(b))
              Z = expecttimes(b,a.*c); 
       else
           Z = expecttimes(a,b.*c); 
       end
      end
  else
      error('pas defini')   
  end        
    
        
end

