function disp(ls,varargin)
% function disp(ls,varargin)

switch class(ls)
    case 'LEVELSET'
        fprintf(['(' class(ls) ' #' num2str(getnumber(ls)) ')'])
    otherwise
        fprintf(['(LEVELSET #' num2str(getnumber(ls)) ' -> ' class(ls) ')'])
end

if ~iseval(ls)
    
    if isa(ls.value{1,1},'function_handle')
        levels = ls.value{1,2};
        disp(' ')
        fprintf('    function_handle : %s(ls1,...)',func2str(ls.value{1,1}))
        disp(' ')
        for k=1:length(levels)
            fprintf('ls%d = ',k);
            disp(levels{k},'novalue','noenrich')
        end
        
    else
        for k=1:size(ls.value,1)
            if israndom(ls.value{k,2})
                ls.value{k,2} = num2str(ls.value{k,2});
            end
        end
        disp(' ')
        ls.value = ls.value';
        disp(struct(ls.value{:}))
        
    end
    
    if ~ischarin('noenrich',varargin)
        disp(struct('enrich',ls.enrich,'enrichtype',ls.enrichtype))
    end
    
    
else
    
    disp(' ')
    if nargin==1 || ~ischarin('novalue',varargin)
        disp(ls.value)
        if ~ischarin('noenrich',varargin)
            disp(struct('enrich',ls.enrich,'enrichtype',ls.enrichtype));
        end
    else
        
        if ~ischarin('noenrich',varargin)
            disp(struct('value',ls.value,'enrich',ls.enrich,'enrichtype',ls.enrichtype))
        else
            disp(struct('value',ls.value));
        end
    end
end

if ~isempty(ls.material)
    disp(ls.material)
end


