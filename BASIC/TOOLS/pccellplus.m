function [U,L]=pccellplus(V,M,k)
% function [u,L]=pccellplus(v,M,k)
% v cell contenant des MYDOUBLEND
% M cell contenant des PCMATRIX ou des POLYCHAOS
% k entier indiquant la dimension du MYDOUBLEND correspondant au
% stochastique
% si M est un POLYCHAOS, la dimension k de v correspond aux coefficients
% sur le chaos
% si M est une PCMATRIX, la dimension k de v correspond aux coefficients
% sur la base de variables al�atoires que sont les composantes de la PCMATRIX
[rep1,pos1]=isclassin(M,'PCMATRIX');
[rep2,pos2]=isclassin(M,'POLYCHAOS');




if rep1 & rep2
    for i=pos2
        s = size(V{i});
        n = length(s)+1;
        v=permute(double(M{i}),[]);


        if isclassin(M,'POLYCHAOS') & isclassin(M,'PCMATRIX')

            for k=1:length(L)

                u = v{1}
                L = M{1};



