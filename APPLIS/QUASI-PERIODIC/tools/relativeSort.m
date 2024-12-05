function [sorted,match] = relativeSort(unsorted,relative,reference,tol,sizeAsRef)
% sorted = relativeSort(unsorted,relative,reference,tol,sizeAsRef)

if nargin < 5
    sizeAsRef = false ;
    if nargin < 4
        tol = 1e-9 ; % tolerance on relative and reference
    end
end

%% Cell recursion
if iscell(unsorted)
    if isnumeric(relative)
        relative = repmat({relative},size(unsorted)) ;
    end
    if isnumeric(reference)
        reference = repmat({reference},size(unsorted)) ;
    end
    sorted = cell(size(unsorted)) ;
    for i = 1:numel(unsorted)
        sorted{i} = relativeSort(unsorted{i},relative{i},reference{i},tol) ;
    end
    return
end

% Safety
if isempty(unsorted)
    sorted = [] ;
    return
end

%% Method

% Round to tolerance
roundedRel = round(relative/tol)*tol ;
roundedRef = round(reference/tol)*tol ;
[match,transfer] = ismember(roundedRef,roundedRel,'rows') ;
if ~sizeAsRef && ~all(match)
    warning('relative does not match perfectly reference.') ;
end
transfer = transfer(match) ;

szu = size(unsorted) ;
if szu(1) == szu(2) % operator
    sorted = unsorted(transfer,transfer) ;
elseif szu(1) == length(relative)
    sorted = unsorted(transfer,:) ;
else
    sorted = unsorted(:,transfer) ;
end

if sizeAsRef
   shortSorted = sorted ;
   sorted = zeros(size(reference,1),size(shortSorted,2)) ;
   sorted(match) = shortSorted ;
end

end
%TODO: implement interpolation across meshes ?