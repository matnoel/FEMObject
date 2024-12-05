function u = MULTIMYDOUBLEND(a,varargin)
% function u = MULTIMYDOUBLEND(a,multidim,sm)
% a : MYDOUBLEND
% multidim : indique la dimension multi dans a
% sm : taille au niveau multi

if nargin==0

    u.multidim = 5;
    u.value = zerosND(1,1,1,1,0);     

    u=class(u,'MULTIMYDOUBLEND');
    superiorto('MYDOUBLEND');

elseif nargin==2 

    u.multidim = varargin{1};

    if length(u.multidim)>1
        error('la dimension stochastique doit etre un entier')
    end

    if isa(a,'MULTIMATRIX')
        u.sm = a.sm;
        a = MYDOUBLEND(a,u.multidim);
    elseif isa(a,'cell')
        [rep,pos] = isclassin('MYDOUBLEND',a(:));
        if length(pos)~=numel(a)
            error('les cellules doivent etre des MYDOUBLEND')
        end
        if nargin>=3   
            u.sm = varargin{2};
        else
            u.sm = size(a);
        end    

    elseif ~isa(a,'MYDOUBLEND')  
        error('rentrer un MYDOUBLEND ou une MULTIMATRIX')
    else
        if nargin>=3   
            u.sm = varargin{2};
        else
            u.sm = [size(a,u.multidim),1];
        end

    end    

    if isa(a,'MYDOUBLEND')
        s = size(a);
        if u.multidim>length(s)
            s=[s,ones(1,u.multidim-length(s))];
        end
        if prod(u.sm)~=s(u.multidim)
            error(['le dimension ' num2str(u.multidim) ' du MYDOUBLEND doit correspondre � celle de ' class(varargin{1})])
        end
    end

    u.value = a;
    u=class(u,'MULTIMYDOUBLEND');
    superiorto('MYDOUBLEND');

else
    error('mauvais argument')
end

