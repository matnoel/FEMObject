function u = PCTIMEMATRIX(v,T,s)
% function u = PCTIMEMATRIX(v,T,s)
% v : PCMATRIX ou PCRADIALMATRIX
% T : TIMEMODEL
% s : taille de la matrice
%
% La deuxieme dimension correspond aux composantes
% temporelles

if nargin>=3 && length(s)~=2
    error('rentrer 2 composantes ou la taille de la matrice');
end


if nargin==0
    u.value=[];
    u.s = size(u.value);
    T = TIMEMODEL();     
    u=class(u,'PCTIMEMATRIX',T);
    superiorto('TIMEMATRIX','PCRADIALMATRIX','PCMATRIX','MULTIMATRIX');

elseif isa(v,'PCTIMEMATRIX')
    if nargin>=3
        u=reshape(v,s);
    else
        u=v;
    end

else


    T = gettimemodel(T);
    T = resetevolparam(T);

    n = length(gettapprox(T));

    if ~israndom(v)
        error('l''argument n''est pas aleatoire, utiliser TIMEMATRIX')  
    end

    if ~israndom(v)
        error('l''argument n''est pas aleatoire')
    end

    if isa(v,'cell') 
        if length(v)~=n
            error('doit avoir le bon nombre de cellules ')
        end

%if ~isclassin('PCRADIALMATRIX',v)
%try
%s = size(v{1});
%u = TIMEMATRIX([v{:}],T,s);
%return 
%end
        s = size(v{1});  
        u.value = v;
        u.s = s;

    else 

        if nargin<3 
            if size(v,2)~=n
                error('la dimension 2 doit correspondre avec l''aproximation temporelle')
            end
            s = size(v);
            s(2)=1;
        end


        if n*prod(s)~=prod(size(v)) || size(v,2)~=n
            error('la matrice aleatoire n''a pas le bon nombre de composantes')
        end

        u.value = v;  
        u.s = s;

    end
    u=class(u,'PCTIMEMATRIX',T);
    superiorto('TIMEMATRIX','PCRADIALMATRIX','PCMATRIX','MULTIMATRIX');


end

