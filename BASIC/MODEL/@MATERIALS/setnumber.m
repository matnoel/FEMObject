function u = setnumber(u,k,l)
% function u = setnumber(u,k)
% renumerotation des materiaux
% signale une erreur si deux materiaux ont le meme numero
% k de taille u.n
%
% function u = setnumber(u,k,l)
% renumerotation des materiaux k
% l sont les nouveaux numeros
%
% function u = setnumber(u)
% numerotation des materiaux sans numero


if nargin==1
    number = getnumber(u);
    
    for k=1:u.n
        if isempty(number{k})
            number{k} = newnumber(u);
            u.MAT{k} = setnumber(u.MAT{k},number{k});
        end
    end
    
elseif nargin==2
    if isa(k,'double')
        k = num2cell(k);
    end
    
    if length(k)~=u.n
        error('All the materials must be renumbered.')
    end
    for j=1:u.n
        u.MAT{j} = setnumber(u.MAT{j},k(j));
    end
    
elseif nargin==3
    if length(k)~=length(l)
        error('k and l must have the same size.')
    end
    for j=1:length(k)
        u.MAT{k(j)} = setnumber(u.MAT{k(j)},l(j));
    end
    
else
    error('Wrong input arguments')
end
