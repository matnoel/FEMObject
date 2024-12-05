function u = setnumber(u,k,l)
% function u = setnumber(u,k)
% renumerotation des polynomes
% si plusieurs polynomes ont le meme numero, on ne garde que la derniere
% k de taille u.M
%
% function u = setnumber(u,k,l)
% renumerotation des polynomes k
% l sont les nouveaux numeros
%
% function u = setnumber(u)
% numerotation des polynomes sans numero

if nargin==1
    number = getnumber(u);
    
    for k=1:u.M
        if isempty(number{k})
            number{k} = newnumber(number);
            u.h{k} = setnumber(u.h{k},number{k});
        end
    end
    
elseif nargin==2
    if isa(k,'double')
        k = num2cell(k);
    end
    if length(k)~=u.M
        error('tous les polynomes doivent etre renumerotes')
    end
    for j=1:u.M
        u.h{j} = setnumber(u.h{j},k{j});
    end
    
elseif nargin==3
    if length(k)~=length(l)
        error('k et l doivent etre de meme taille')
    end
    for j=1:length(k)
        u.h{k(j)} = setnumber(u.h{k(j)},l(j));
    end
    
else
    error('mauvais arguments')
end


function n = newnumber(number)
% function n = newnumber(number)

allnumber = [number{:}];
n = min(setdiff(1:length(number)+1,allnumber));

return