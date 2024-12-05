function u = setnumber(u,k,l)
% function u = setnumber(u,k)
% renumerotation des variables
% si plusieurs variables ont le meme numero, on ne garde que la derniere
% k de taille u.M
%
% function u = setnumber(u,k,l)
% renumerotation des variables k
% l sont les nouveaux numeros
%
% function u = setnumber(u)
% numerotation des variables sans numero

if nargin==1
    number = getnumber(u);
    
    for k=1:u.M
        if isempty(number{k})
            if isa(u.RV{k},'CONDRANDVAR')
                number{k} = newnumber(u,getstodim(u.RV{k}));
            else
                number{k} = newnumber(u);
            end
            u.RV{k} = setnumber(u.RV{k},number{k});
        end
    end
    
elseif nargin==2
    if isa(k,'double')
        k = num2cell(k);
    end
    
    if length(k)~=u.M
        error('toutes les variables doivent etre renumerotees')
    end
    for j=1:u.M
        u.RV{j} = setnumber(u.RV{j},k{j});
    end
    
elseif nargin==3
    if isa(k,'double')
        k = num2cell(k);
    end
    if isa(l,'double')
        l = num2cell(l);
    end
    
    if length(k)~=length(l)
        error('k et l doivent etre de meme taille')
    end
    for j=1:length(k)
        u.RV{k{j}} = setnumber(u.RV{k{j}},l{j});
    end
    
else
    error('mauvais arguments')
end

