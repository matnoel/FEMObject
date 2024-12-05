function param = getparam(mat,i)
% function param = getparam(mat,i)

if nargin==1
    param = mat.param;
else
    if isa(i,'double')
        param = mat.param(i,2);
    elseif isa(i,'char') || isa(i,'cell')
        [rep,pos] = ischarin(i,mat.param(:,1));
        param = mat.param(pos(find(rep)),2);
    end
end

if length(param)==1
    param = param{1};
elseif isempty(param)
    param = [];
end

