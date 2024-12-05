function ls=lseval(ls,D,varargin)

tolls = getfemobjectoptions('tolerancelevelset');
if iseval(ls)
    return
elseif nargin==1
    error('rentrer un MODEL ou un POINT en arguments en argument')
elseif israndom(ls)
    
    if nargin<3 || ~isclassin('POLYCHAOS',varargin)
        error('rentrer un POLYCHAOS en 3e argument')
    end
    PC = getPC(getclassin('POLYCHAOS',varargin));
    
    
    value = decompmatrix(PC,[],[],@(a) getvalue(lseval(randomeval(ls,a,RANDVARS(PC)),D)));
    
    ls = setvalue(ls,value);
    
elseif isa(ls.value{1,1},'function_handle')
    
    fun = fcnchk(ls.value{1,1});
    levels = ls.value{1,2};
    
    for k=1:length(levels)
        levelsD{k}=lseval(levels{k},D,varargin{:});
    end
    lsini = ls;
    ls = fun(levelsD{:});
    ls.number=lsini.number;
    ls.material = lsini.material;
    ls.enrich = lsini.enrich;
    ls.enrichtype = lsini.enrichtype;
    
elseif isa(D,'MODEL')
    
    P = POINT(getnode(D));
    
    lsx = ls.sign*double(eval(ls,P,ls.value{:,2}));
    lsx(find(abs(lsx)<tolls))=0;
    value = MULTIMATRIX(lsx(:),[getnbnode(D),1]);
    
    if isa(value,'MULTIMATRIX') && length(value)==1 & ~israndom(value)
        value = double(value);
    end
    
    ls = setvalue(ls,value);
    ls.sign=1;
    
elseif isa(D,'NODE')
    lsx = lseval(ls,POINT(D));
    lsx = double(lsx);
    lsx = reshape(lsx,getnbnode(D),size(lsx,4));
    if ischarin('LEVELSET',varargin)
        value = MULTIMATRIX(lsx(:),[getnbnode(D),1]);
        if isa(value,'MULTIMATRIX') && length(value)==1 && ~israndom(value)
            value = double(value);
        end
        ls = setvalue(ls,value);
        ls.sign=1;
    else
        ls=lsx;
    end
else
    P = POINT(D);
    ls = ls.sign*eval(ls,P,ls.value{:,2});
    if tolls>0
        ls = double(ls);
        ls(find(abs(ls)<tolls))=0;
        ls = MYDOUBLEND(ls);
    end
end

if isa(ls,'LEVELSET') && ~israndom(ls) && isclassin('POLYCHAOS',varargin)
    %PC = getPC(getclassin('POLYCHAOS',varargin));
    %ls.value = full(expand(ls.value*one(PC)));
end

function value = funeval(a,ls,P,tolls)
ls.value = funrandomparam(ls.value,@randomeval,a);

lsx = ls.sign*double(eval(ls,P,ls.value{:,2}));
lsx = sparse(lsx(:));
lsx(find(abs(lsx)<tolls))=0;
value = MULTIMATRIX(lsx(:),[numel(P),1]);

return