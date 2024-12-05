function rvs = RANDVARS(varargin)
% function rvs = RANDVARS(varargin)
% ATTENTION : toutes les variables aleatoires sont considerees independantes
% si plusieurs variables ont le meme numero, on ne garde que la derniere

if nargin==1 && isa(varargin{1},'RANDVARS')
    rvs = varargin{1};
elseif nargin==2 && isa(varargin{1},'RANDVAR') && isa(varargin{2},'double')
    r=cell(1,varargin{2});
    for i=1:varargin{2}
        r{i} = varargin{1};
    end
    rvs = RANDVARS(r{:});
    
else
    rvs.RV=cell(1,0);
    for k=1:nargin
        if isa(varargin{k},'RANDVAR') || isa(varargin{k},'CONDRANDVAR')
            rvs.RV = [rvs.RV , {varargin{k}}];
        elseif isa(varargin{k},'RANDVARS')
            rvs.RV = [rvs.RV , varargin{k}.RV];
        elseif isa(varargin{k},'cell')
            rvsk = RANDVARS(varargin{k}{:}) ;
            rvs.RV = [rvs.RV , rvsk.RV];
        end
    end
    rvs.M = length(rvs.RV);
    rvs = class(rvs,'RANDVARS');
    rvs = setnumber(rvs);
    for k=1:rvs.M
        if isa(rvs.RV{k},'CONDRANDVAR')
            
            Y = getY(rvs.RV{k});
            stodim = getstodim(rvs.RV{k});
            for j=1:length(Y)
                if isempty(stodim{j}) || ~ismember(Y{j},rvs)
                    if isempty(stodim{j})
                        stodim{j} = newnumber(rvs);
                        Y{j} = setnumber(Y{j},stodim{j});
                        rvs.RV{k} = setstodim(rvs.RV{k},stodim);
                    end
                    rvs = RANDVARS(rvs,Y{j});
                else
                    [ok,rep] = ismember(Y{j},rvs);
                    
                    if Y{j}~=rvs.RV{rep}
                        error('deux variables de meme numero sont differentes')
                    end
                end
            end
        end
    end
    rvs = unique(rvs);
end
