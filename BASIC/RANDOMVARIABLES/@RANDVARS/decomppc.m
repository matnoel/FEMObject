function [apc,H,p]=decomppc(a,varargin)
% function [apc,H,p]=decomppc(a,'order',p,'tol',tol)
% decomposition de a sur le chaos polynomial
% p : ordre de la decomposition
% tol : erreur souhaitee
%
% function apc=decomppc(a,'order',p,H)
% H : RANDPOLYS -> decomposition sur le chaos generalise ou element finis
% H{i} est un RANDPOLY pour la dimension i (eventuellement un POLYFE)
%
% function apc=decomppc(a,'fedim',i,'femesh',n)
% FE au niveau stochastique pour les dimensions i 
% n : nombre de subdivisions de [0,1]
% n : tableau de (length(i)) cellules. n{j} contient le decoupage de [0,1]
% pour la dimension i(j)
%
% function apc = PCMODEL(a,'order',p,'pcgdim',i)
% PC generalise au niveau stochastique pour les dimensions i 
% appel de RANDPOLY(a{i}) pour determiner la base polynomiale de la VA a{i}


p = getcharin('order',varargin);
if ~isempty(p) & length(p)==1
p=repmat(p,1,a.M);
end
tol=getcharin('tol',varargin);
if ~isempty(tol) & length(tol)==1
tol=repmat(tol,1,a.M);
end

H = getclassin('RANDPOLYS',varargin);

if ~isempty(H)
    if length(H)~=a.M
error('RANDPOLYS doit contenir autant de polynomes que le nombre de VA')
    end
else

H=RANDPOLYS();

H(1:length(a))=POLYHERMITE();
pcgdim = getcharin('pcgdim',varargin);

H(pcgdim) = RANDPOLYS(a.RV(pcgdim));
fedim = getcharin('fedim',varargin);
if ~isempty(fedim)
n=getcharin('femesh',varargin); 
if ~isa(n,'cell') || length(n)~=length(fedim)
error('femesh must be a cell array')
end
for k=1:length(fedim)
H(fedim(k)) = POLYFE(n{k},p(k));
end
end

end

H = setnumber(H,getnumber(a));
apc = cell(1,a.M);
for k=1:a.M

   options={};
   if ~isempty(p)
   options=setcharin('order',options,p(k));
   end
   if ~isempty(tol)
   options=setcharin('tol',options,tol(k));
   end

   if isa(H{k},'POLYFE')
   apc{k} = decompfe(a.RV{k},H{k},options{:});
   else
   apc{k} = decomppc(a.RV{k},H{k},options{:});         
   end

end
