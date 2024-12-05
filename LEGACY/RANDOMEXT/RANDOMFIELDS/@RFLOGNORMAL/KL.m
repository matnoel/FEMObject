function [RFpc,varargout]=KL(RF,m,S,varargin)
% function RFpc=KL(RF,m,S)
% RF : RFLOGNORMAL
% S : MODEL de type MAILLAGE
% m : ordre de la decomposition
% RFpc : DISCRANDFIELD 
%
% function RFpc=KL(RF,m,S,'pcorder',p)
% p : ordre du chaos polynomial pour la decomposition des variables
% aleatoires
%
disp('---- Decomposition du champ lognormal----')

S=createddlnode(S,DDL('U'));
x = getcoord(S.node);
marg = marginal(RF,x);
mu = getparam(marg,'mu');
sigma = getparam(marg,'sigma');
Gam = RFGAUSSIAN(mu,sigma,RF.correl);

V = KL(Gam,m,S,'rescale','modes');

V0 = exp(mu+sigma.^2/2);
V=V(:,2:end);


p = getcharin('pcorder',varargin,1);

fprintf(' -> decomposition sur le chaos : degre p=%d\n' , p )
PC = POLYCHAOS(m,p);


n = size(V,1);
RFpc = zeros(n,length(PC));
indices = getindices(PC);
for k=1:size(indices,1)
   pourcentage(k,size(indices,1),100);
   ind = indices(k,1:end-1); 
   fact = prod(sqrt(factorial(ind))); 
   prodV = prod(V.^(repmat(ind,n,1)),2);
   RFpc(:,k) = V0.*prodV/fact;
end


RFpc = PCMATRIX(RFpc,[n,1],PC);



if nargout>1 && ischarin('modes',varargin)
    varargout{1}=V;
end