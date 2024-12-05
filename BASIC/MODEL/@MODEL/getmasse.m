function masse=getmasse(M,derivation)

% Masse :
if nargin==1
    derivation=[0 0 0];
end
masse=setfree(MULTILINFORM(derivation,1,0),0);
masse=masse{M}(:,:,:);