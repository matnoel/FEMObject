function x = MULTIMATRIX(a,s,sm)
% function x = MULTIMATRIX(a,s,sm)
% s : taille d'une matrice
% a : double de taille prod(s)*(nombre de matrices)

% function x = MULTIMATRIX(a,k)
% a : MYDOUBLEND
% k : dimensions correspondant aux MULTIMATRIX


if nargin==0
    x.s=[0,0];
    x.value=sparse(0,0);
    x.sm = [0,1];

    x=class(x,'MULTIMATRIX');


elseif isa(a,'MULTIMATRIX') && nargin==1
    x = a;

elseif isa(a,'MULTIMATRIX') 
    issames = (length(s)==length(a.s) && all(s==a.s));
    if nargin==3
        issamesm = (length(sm)==length(a.sm) && all(sm==a.sm));
    end
    if (nargin==2 && issames) || (nargin==3 && issames && issamesm)        
        x = a;
    elseif nargin==3 && prod([a.s,a.sm])==prod([s,sm])
        if iscell(a)
            x=a;
            if prod(s)<prod(a.s)
                a.value = reshape(a.value,1,prod(a.sm));
                x.value = cell(sm);
                for k=1:length(a.value)
                    a.value{k} =  mat2cell(a.value{k},prod(s),ones(1,prod(sm)/prod(a.sm)));   
                    x.value(:,k) = a.value{k};
                end
                x = reshapem(x,sm);
                x=reshape(x,s);

            else
                error('pas programme')    
            end

        else     
            x = MULTIMATRIX(double(a.value),s,sm);  

        end

    else
        error('pas prevu')
    end

elseif isa(a,'MULTIMATRIX') && nargin==3 && length(s)==length(a.s) && all(s==a.s)

elseif isa(a,'cell')
    okcell=1;

    if nargin<2    
        s=size(a{1});
    else
        for k=1:length(a)
            a{k}=reshape(a{k},s);
        end
    end

    if nargin<3
        sm = size(a);
    else
        a = reshape(a,sm);
    end
    x.s=s;
    if okcell
        x.value=a;
    else
        x.value = reshape([a{:}],prod(s),prod(sm));    
    end
    x.sm = sm;

    x=class(x,'MULTIMATRIX');

elseif isa(a,'MYDOUBLEND')
    k=s;
    s=size(a);
    s=[s,ones(1,max(k)-length(s))];
    nk = setdiff(1:length(s),k);
    sm = [s(k),1];
    s = s(nk);
    a = permute(double(a),[nk,k]);
    x = MULTIMATRIX(a,s,sm);

elseif isa(a,'double') 

    if nargin<2  
        s = size(a);
        sm=[1,1];
    elseif nargin<3
        if prod(s)==0
            error('MULTIMATRIX indeterminee car 0 comoposantes dans la matrice') 
        end
        sm = [numel(a)/prod(s),1];
    end

    x.s = s ;
    x.value = reshape(a,[prod(s) prod(sm)]);
    x.sm = sm;

    x=class(x,'MULTIMATRIX');

else
    error('')
end



if numel(x.sm)==1
    error('un seul element dans la dimension multi de la multimatrix')
end
