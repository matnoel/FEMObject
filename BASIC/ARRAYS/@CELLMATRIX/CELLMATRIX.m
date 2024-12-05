function x = CELLMATRIX(a,s,sm)
% function x = CELLMATRIX(a,s,sm)
% s taille d'une matrice
% a : double de taille prod(s)*(nombre de matrices)

% function x = CELLMATRIX(a,k)
% a : MYDOUBLEND
% k : dimensions correspondant aux CELLMATRIX


if nargin==0
    x.s=[0,0];
    x.value=cell(0,0);
    x.sm = 0;

    x=class(x,'CELLMATRIX');


elseif isa(a,'CELLMATRIX')
    x = a;

elseif isa(a,'cell')

    if nargin<2 | isempty(s)    
        s=size(a{1});
%for k=1:numel(a)
%    a{k}=reshape(a{k},s);
%end
    end
    if nargin<3 | isempty(sm)
        a = reshape(a,sm);
    end
    x.s = s;
    x.value = a ;
    x.sm = sm ;

    x=class(x,'CELLMATRIX');

elseif isa(a,'MYDOUBLEND')
    k=s;
    s=size(a);
    s=[s,ones(1,max(k)-length(s))];
    nk = setdiff(1:length(s),k);
    sm = s(k);
    s = s(nk);
    a = permute(double(a),[nk,k]);
    x = CELLMATRIX(a,s,sm);

elseif isa(a,'double') & nargin==1
    x.s = size(a);
    x.value = {a};
    x.sm = [1,1];
    x=class(x,'CELLMATRIX');

else

    x.s = s;
    if nargin<3
        if prod(s)==0
            error('CELLMATRIX indeterminee car 0 comoposantes dans la matrice') 
        end
        sm = [numel(a)/prod(s),1];
    else
        if length(sm)==1
            sm=[sm,1];
        end
    end

    x.value = reshape(a,[prod(s) prod(sm)]);
    m = prod(sm);
    n = prod(s);
    x.value = mat2cell(a,[repmat(n,m,1),repmat(m,n,1)]);

    x.sm=sm;


    x=class(x,'CELLMATRIX');

end

