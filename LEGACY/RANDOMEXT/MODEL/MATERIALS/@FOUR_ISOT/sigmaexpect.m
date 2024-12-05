function [se,B] = sigmaexpect(mat,elem,xnode,xgauss,qe,l1,l2,varargin)
% function [se,B] = sigmaexpect(mat,elem,xnode,xgauss,qe,l1,l2,varargin)

B = calc_B(elem,xnode,xgauss);
gradu = B*qe;

k = evalparampc(mat,'k',[],elem,xnode,xgauss);
se = k*gradu;

if ~israndom(se)
    se = se*expect(l1,l2);
else
    se = expect(se,l1,l2);
end

if mat.b || mat.k2 || mat.r
    N = calc_N(elem,xnode,xgauss);
    u = N*qe;
    
    if mat.b
        b = evalparampc(mat,'b',[],elem,xnode,xgauss);
        seb = b*u;
        if ~israndom(seb)
            seb = seb*expect(l1,l2);
        else
            seb = expect(seb,l1,l2);
        end
        se = se + seb;
    end
    
    if mat.k2
        k2 = evalparampc(mat,'k2',[],elem,xnode,xgauss);
        sek2 = k2*(u.*u)*gradu;
        if ~israndom(sek2)
            sek2 = sek2*expect(l1,l2);
        else
            sek2 = expect(sek2,l1,l2);
        end
        se = se + sek2;
    end
    
    if mat.r
        r = evalparampc(mat,'r',[],elem,xnode,xgauss);
        ser = r*(u.*u)*u;
        if ~israndom(ser)
            ser = ser*expect(l1,l2);
        else
            ser = expect(sea,l1,l2);
        end
        se = se + ser;
    end
end


