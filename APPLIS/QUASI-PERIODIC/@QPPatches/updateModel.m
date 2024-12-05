function ps = updateModel(ps,model)
% ps = updateModel(ps,model)

if nargin == 1
    model = getModel(ps) ;
end

cA = updateModel(getConductivityAssembler(ps),model) ;
ps = updateConductivityAssembler(ps,cA) ;

end

