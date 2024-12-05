function X = randomknowing(rv,y,RV,varargin)
% function X = randomknowing(rv,y,RV,varargin)
if ~isa(y,'cell')
    error('rentrer une cellule pour y')
end
    
alleval = 1;
for i=1:length(rv.Y)
if isa(rv.Y{i},'RANDVAR') || isa(rv.Y{i},'CONDRANDVAR')
[rep,loc]= ismember(rv.Y{i},RV) ;  

if rep
rv.Y{i} = y{loc};
else
alleval = 0;
end
end
end

if alleval
for i=1:length(rv.funparam)
    if isa(rv.funparam{i},'inline') || isa(rv.funparam{i},'function_handle')
    rv.funparam{i}=rv.funparam{i}(rv.Y{:});
    end
end
    X = setnumber(rv.X(rv.funparam{:}),rv.number); 
    
    X = random(X,varargin{:});
else
    X=rv;
end
