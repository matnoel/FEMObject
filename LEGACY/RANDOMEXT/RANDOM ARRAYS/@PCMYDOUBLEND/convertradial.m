function a = convertradial(a)
% function u = convertradial(a)
% a : PCMYDOUBLEND

if ~isradial(a)

    if numel(a)==1
   
    L = double(a.V);
    L = PCMATRIX(L(:),[1,1],getPC(a));
   
    a = PCMYDOUBLEND(MYDOUBLEND(1),L,a.stodim);
    else
     
      L = MULTIMATRIX(a.V,a.stodim);
      L = double(L);
      n = size(L,1);
      L = PCMATRIX(L,[n,1],getPC(a));
      
      V = diag(ones(1,n));
      s = [size(a),n];
      V = reshape(V,s);
      p = 1:max(length(s),a.stodim);
      p(a.stodim) = length(s);
      p(length(s))=a.stodim;
      V = MYDOUBLEND(permute(V,p));
      a = PCMYDOUBLEND(V,L,a.stodim);
    end
    
      
end
