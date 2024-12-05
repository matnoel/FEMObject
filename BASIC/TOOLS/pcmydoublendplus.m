function [U,L]=pcmydoublendplus(V,M,k)
% function [u,L]=pcmydoublendplus(v,M,k)
% v cell contenant des MYDOUBLEND
% M cell contenant des PCMATRIX ou des POLYCHAOS
% k entier indiquant la dimension du MYDOUBLEND correspondant au
% stochastique
% si M est un POLYCHAOS, la dimension k de v correspond aux coefficients
% sur le chaos
% si M est une PCMATRIX, la dimension k de v correspond aux coefficients
% sur la base de variables al�atoires que sont les composantes de la PCMATRIX

for i=1:length(V)
    if isa(M{i},'PCMATRIX')
        V{i} = PCRADIALMATRIX(MULTIMATRIX(V{i},k),size(V{i}),M{i});
    elseif strcmp(class(M{i}),'POLYCHAOS')
        V{i} = PCMATRIX(MULTIMATRIX(V{i},k),M{i});
    elseif isempty(M{i})
        V{i} = double(V{i});    
    else
        error('pas prevu')
    end
    if i==1
        U = V{1};
    else
        U = U + V{i}; 
    end
end

s = size(U);
n = [1:length(s),k];
n(k) = length(s)+1;

if isa(U,'PCRADIALMATRIX')
    L = getD(U)*getL(U);    
    U = permute(MYDOUBLEND(getV(U)),n);
elseif isa(U,'PCMATRIX')
    L = getPC(U);    
    U = permute(MYDOUBLEND(getmultimatrix(U)),n);    
elseif isa(U,'double')
    L = []; 
    U = MYDOUBLEND(U);
end

