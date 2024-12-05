function b0 = create_initial_condition(T,g0)
% function b0 = create_initial_condition(T,g0)

if ~isa(T,'DGTIMESOLVER')
    error('ca marche avec DG')
end
    
f0 = getinitialequivalent(T);
b0 = g0*f0;
