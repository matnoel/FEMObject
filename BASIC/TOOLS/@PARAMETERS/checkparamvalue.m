function varargout = checkparamvalue(P,param,varargin)
% function P = checkparamvalue(P,'paramname',paramvalue1,paramvalue2)
% verifient que le parametre de nom param vaut paramvalue1 ou ...
% si ce n'est pas le cas, message d'erreur et affichage des possibilites

paramvalues = varargin;
paramvalue = getfield(P.param,param);
rep = 0 ;

switch class(paramvalue)
    case 'char'
        for i=1:length(paramvalues)
            if isa(paramvalues{i},'char') && strcmp(paramvalues{i},paramvalue)
                rep=1;
                break
            end
        end
    case {'double','logical'}
        if ndims(paramvalue)>2
            error('ne verifie pas les ND-array');
        end
        for i=1:length(paramvalues)
            if (isa(paramvalues{i},'double') || isa(paramvalues{i},'logical') ) && ...
                    all(size(paramvalues{i})==size(paramvalue)) && ...
                    all(all(paramvalues{i}==paramvalue))
                rep=1;
                break
            end
        end
    otherwise
        error('ne verifie pas ce type d''objets')
        
        for i=1:length(paramnames)
            if ~isfield(P.param,paramnames{i})
                P.param = setfield(P.param,paramnames{i},paramvalues{i});
            end
        end
        
end

if rep==0
    fprintf('Mauvaise valeur du parametre "%s". Les possiblites sont : \n',param);
    disp(varargin)
    error(' ')
end

if nargout==1
    varargout{1}=rep;
end

