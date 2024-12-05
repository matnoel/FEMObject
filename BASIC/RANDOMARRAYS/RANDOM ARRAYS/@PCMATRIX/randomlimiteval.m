function a=randomlimiteval(apc,x,RV,varargin)
% function a = randomeval(apc,x)
% evaluation de la PCMATRIX pour les realisations x
% M variables aleatoires
% x : double n-by-M (n : nombre de realisations) 
% ou x : 1-by-M cell de n-by-1 double 
%
% a : MULTIMATRIX de taille size(apc) et de longueur n si n>1 et numel(apc)>1
% a : double sinon
%
% function a = randomeval(apc,x,RV)
% RV indique les dimensions stochastiques de x (pour la correspondance avec apc)
% -> utilisation de [ok,rep] = ismember(apc,RV) et x = x(:,rep)
% RV peut etre un RANDPOLYS ou une RANDVARS ou encore un double
if nargin == 2
    RV = RANDVARS(apc.POLYCHAOS);
    if (isa(x,'cell') & length(x)~=length(RV)) || (isa(x,'double') & size(x,2)~=length(RV))
     error(['on attend ' num2str(length(RV)) ' jeux de valeurs' ])
    end
    if isa(x,'cell')
    x=[x{:}];    
    end
end

if nargin>=3
x = transfer(RANDVARS(RV),RANDVARS(apc.POLYCHAOS),x);
end 
%keyboard
  nmax = getcharin('nmax',varargin,10);
  n = ceil(size(x,1)/nmax);
  nb = floor(size(x,1)/n);
  a = [];
  s = 0;
  for i=1:nb
      H=polyval(apc.POLYCHAOS,x((i-1)*n+1:i*n,:));
      s = s + size(H,1); 
      if iscell(apc)    
      a =[a getvalue(multimtimes(H,apc.MULTIMATRIX))];
      else
      a = [a double(apc) * H']; 
      end    
  end
  if size(x,1)>n*nb
     H=polyval(apc.POLYCHAOS,x(n*nb+1:end,:)); 
     if iscell(apc)    
      a =[a getvalue(multimtimes(H,apc.MULTIMATRIX))];
     else
      a = [a double(apc) * H']; 
      end    
  end
  
  if s==1
  if iscell(apc)  
  a=a{1};  
  end
  a=reshape(a,size(apc));
  elseif all(size(apc)==1)  
  if iscell(apc)
  a = [a{:}];    
  end  
  a=reshape(a,1,s);    
  else
  a=MULTIMATRIX(a,size(apc),[size(x,1),1]); 
  end
