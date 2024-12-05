function assembler = homogeneous(value,model)
% assembler = homogeneous(value,model)
% Static method to create assembler of homogeneous conductivity.

pattern = struct('name','uniform','value',value,'size',[],'center',[],'offset',[]);
assembler = QPConductivityAssembler('model',model,'patterns',pattern,...
    'distribution',{(1:getCellNb(model))'}) ;

end