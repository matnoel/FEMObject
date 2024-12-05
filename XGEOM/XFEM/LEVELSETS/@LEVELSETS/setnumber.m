function u = setnumber(u,k,l)
% function u = setnumber(u,k)
% renumerotation des levelsets
% signale une erreur si deux levelsets ont le meme numero
% k de taille u.n
%
% function u = setnumber(u,k,l)
% renumerotation des levelsets k
% l sont les nouveaux numeros
%
% function u = setnumber(u)
% numerotation des levelsets sans numero

if nargin==1
    number = getnumber(u);
    
    for k=1:u.n
        if isempty(number{k})
            number{k} = newnumber(u);
            u.LS{k} = setnumber(u.LS{k},number{k});
        end
    end
    
elseif nargin==2
    if isa(k,'double')
        k = num2cell(k);
    end
    
    if length(k)~=u.n
        error('toutes les levelsets doivent etre renumerotees')
    end
    for j=1:u.n
        u.LS{j} = setnumber(u.LS{j},k(j));
    end
    
elseif nargin==3
    if length(k)~=length(l)
        error('k et l doivent etre de meme taille')
    end
    for j=1:length(k)
        u.LS{k(j)} = setnumber(u.LS{k(j)},l(j));
    end
    
else
    error('mauvais arguments')
end


