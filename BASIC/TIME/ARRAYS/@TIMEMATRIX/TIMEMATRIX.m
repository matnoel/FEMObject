function u = TIMEMATRIX(v,T,s,varargin)
% function u = TIMEMATRIX(v,T,s)
% v : 'cell' ou double ou PCMATRIX ou PCRADIALMATRIX
% T : TIMEMODEL
% s : taille de la matrice
% storage : pour faire un changement de type de stockage
%
% La deuxieme dimension correspond aux composantes
% temporelles

switch nargin
case 0
    u.value =[];
    u.s = size(u.value);
    T = TIMEMODEL();
    u=class(u,'TIMEMATRIX',T);
    superiorto('PCRADIALMATRIX','PCMATRIX','MULTIMATRIX');


otherwise
    T = gettimemodel(T);
    %T = resetevolparam(T);
    n = length(gettapprox(T));
    if israndom(v)
        if nargin==2
            u = PCTIMEMATRIX(v,T); 
        elseif nargin>=3
            u = PCTIMEMATRIX(v,T,s); 
        else
            u = PCTIMEMATRIX(v);
        end
        return
    else

        if nargin>=3 && ~isempty(s) && length(s)~=2
            error('rentrer 2 composantes ou la taille de la matrice');
        end

        if isa(v,'cell') 
            if length(v)~=n
                error('doit avoir le bon nombre de cellules ')
            end
            if ischarin('cell',varargin)
                u.value = v;
                u.s = size(v{1});
                u=class(u,'TIMEMATRIX',T);
                superiorto('PCRADIALMATRIX','PCMATRIX','MULTIMATRIX');

            else
                s = size(v{1});   
                u = TIMEMATRIX([v{:}],T,s);
            end

        else
            if isa(v,'TIMEMATRIX')
                if nargin==1
                    u=v;
                elseif nargin>=3    
                    u = reshape(v,s);
                end

                return

            elseif isa(v,'MULTIMATRIX') 
                if size(v,2)~=n
                    error('la dimension 2 doit correspondre avec l''aproximation temporelle')
                end

                u.value = v;
                u.s = size(v);
                u.s(2) = 1;
                if nargin>=2 && ~all(u.s==s)
                    error('la dimension 2 doit correspondre avec l''aproximation temporelle')
                end

            elseif isa(v,'double') 
                if nargin<3
                    if numel(v)~=n
                        s = [size(v,1),1];
                    else
                        s=[1,1];    
                    end
                elseif (prod(s)*n)~=numel(v) || size(v,2)~=n
                    error('la dimension 2 doit correspondre avec l''approximation temporelle')
                end

                u.value = reshape(v,prod(s),n);
                u.s = s;

            else 

                if nargin<3
                    s = [size(v,1),1];
                    if size(v,2)~=n
                        error('les dimensions ne correspondent pas avec l''approximation temporelle')    
                    end
                else
                    v = reshape(v,prod(s),n);    
                end
                u.value=v;
                u.s = s;


            end
            u=class(u,'TIMEMATRIX',T);
            superiorto('PCRADIALMATRIX','PCMATRIX','MULTIMATRIX');

        end

    end
end
