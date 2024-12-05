function varargout = getcharin(s,var,def)
% function var = getcharin(propertyname,var,defaultpropertyvalue)
% var : tableau de cellules contenant des paires du type {'propertyname',propertyvalue}
% propertyname : nom d'une propriete dont on veut la valeur
% defaultpropertyvalue : valeur par defaut si la propriete n'existe pas

if isa(s,'char')
    s={s};
end
if nargin==2
    def=cell(0,length(s));
else
    if ~isa(def,'cell')
        def={def};
    end
end

rep=zeros(1,length(s));pos=zeros(1,length(s));
for j=1:length(s)
    for i=1:length(var)
        if isa(var{i},'char') && strcmpi(var{i},s{j})
            rep(j)=1;
            pos(j)=i;
        end
    end
end

for j=1:length(s)
    if rep(j)
        varargout{j} = var{pos(j)+1};
    elseif nargin==3 && ~isempty(def{j})
        varargout{j} = def{j};
    else
        varargout{j} = [];
        
    end
end
