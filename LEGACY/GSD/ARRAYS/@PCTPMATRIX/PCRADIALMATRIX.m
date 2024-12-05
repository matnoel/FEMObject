function x = PCRADIALMATRIX(x)
% function x = PCRADIALMATRIX(x)

if getnbgroups(x)>1
    error('possible only if one group in POLYCHAOSTP')
end
PC = getpcgroup(x,1);
x = PCRADIALMATRIX(x.phi0,size(x.phi0),PCMATRIX(x.phi{1},[1,1],PC));


