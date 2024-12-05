function x = PCARRAY(a,PC,pcdim)
% function x = PCARRAY(a,PC,pcdim)
% a : array
% PC : POLYCHAOS
% pcdim : indique quelle dimension correspond au chaos polynomial
%         par defaut la derniere

if nargin==1 & isa(a,'PCARRAY')
    x = varargin{1};
elseif isa(a,'PCRADIAL') | isa(a,'PCCELL')
    x=expand(a);    
else
    a=MYDOUBLE(a);
    x.pcdim = ndims(a) ;
    x=class(x,'PCARRAY',a,getPC(PC));
    superiorto('MYDOUBLE','POLYCHAOS');

end

