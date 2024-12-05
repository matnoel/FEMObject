function c = component(model,t,dir)
% c = component(model,t,dir)
% Returns component of t along direction dir as a TuckerLikeTensor.
% If t.space is a TSpaceOperators, dir can have two elements.

% For safety
if t.sz(end) == getNbCellDoF(model)
    warning('Tensor is scalar-valued: I return tensor')
    c = t ;
    return
end

% For operators only
if numel(dir)==2
    % Check whether dir is consistent with t.space
    assert(isa(t.space,'TSpaceOperators'),'Is not an operator')
    % Get requested "column"
    t = t*model.basisVector(dir(2)) ;
end

E = model.basisVector(dir(1)) ;
E = toOperator(E) ;
E.space.spaces{end} = cellfun(@(m)m(dir(1):2:end,:),...
    E.space.spaces{end},'UniformOutPut',false) ;
E = updateAllProperties(E) ;

c = E*t ;

%TODO: Use evalAtIndices along dimension model.order ?
end