function metric=getmetric(M,derivation)

% Metric :
if nargin==1
    derivation=[0 0];
end
metric=setfree(MULTILINFORM(derivation,1,0),0);
metric=metric{M}(:,:);