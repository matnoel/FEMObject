function b = unfreevector(S,b)
% function b = unfreevector(S,b)

if size(b,1)==getnbpoints(S)
    return
end
if isa(b,'TIMEMATRIX')
    btemp = unfreevector(S,double(b));
    b = TIMEMATRIX(btemp,gettimemodel(b),[size(btemp,1),1]);
elseif isa(b,'TIMERADIALMATRIX')
    b = cell2mat(b);
    btemp = unfreevector(S,double(getV(b)));
    b = TIMERADIALMATRIX(btemp,[size(btemp,1),1],getL(b));
    
else
    
    
    if getnbddl(S)~=size(b,1)
        switch S.bc.type
            case 'dirichlet'
                rep1 = S.bc.node;
                rep2 = setdiff(1:getnbddl(S),rep1);
                btemp = b;
                b = zeros(getnbddl(S),size(b,2));
                b(rep2,:) = btemp;
                
            case 'periodic'
                b = [b;b(1,:)];
        end
        
    end
    
end