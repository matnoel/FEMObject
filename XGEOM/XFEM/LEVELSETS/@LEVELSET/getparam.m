function param = getparam(ls,i)
% function param = getparam(ls,i)

if isalevelset(ls)
    error('pas de paramatres pour une LEVELSET')
else
    if nargin==1
        param = ls.value;
    else
        if isa(i,'double')
            param = ls.value(i,2);
        elseif isa(i,'char') || isa(i,'cell')
            [rep,pos] = ischarin(i,ls.value(:,1));
            param = ls.value(pos(find(rep)),2);
        end
    end
    
    if length(param)==1
        param = param{1};
    elseif isempty(param)
        param = [];
    end
end