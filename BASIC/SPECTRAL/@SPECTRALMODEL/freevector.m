function b = freevector(S,b)
% function b = freevector(S,b)

if isa(b,'TIMEMATRIX')
    btemp = freevector(S,double(b));
    b = TIMEMATRIX(btemp,gettimemodel(b),[size(btemp,1),1]);
elseif isa(b,'TIMERADIALMATRIX')
    btemp = freevector(S,double(getV(b)));
    b = TIMEMATRIX(btemp,[size(btemp,1),1],getL(b));
    
else
    if ~isempty(S.bc)
        switch S.bc.type
            case 'dirichlet'
                rep = S.bc.node;
                b(rep,:)=[];
            case 'periodic'
                b(1,:) = b(1,:)+b(end,:);
                b(end,:)=[];
        end
    end
    
end