function [se,B] = sigma(mat,elem,xnode,xgauss,qe,varargin)
% function [se,B] = sigma(mat,elem,xnode,xgauss,qe,varargin)

k = evalparam(mat,'k',elem,xnode,xgauss);
B = calc_B(elem,xnode,xgauss);
gradu = B*qe;

se = k*gradu;

if mat.b || mat.k2 || mat.r || mat.r2 || mat.r3
    N = calc_N(elem,xnode,xgauss);
    u = N*qe;
    
%     if mat.b
%         warning('mal calcule')
%         b = getparam(mat,'b');
%         if isa(b,'cell')
%             bx = b{1};
%             by = b{2};
%             bx = evalparam(setparam(mat,'b',bx),'b',elem,xnode,xgauss);
%             by = evalparam(setparam(mat,'b',by),'b',elem,xnode,xgauss);
%             b = [bx;by];
%         else
%             b = evalparam(mat,'b',elem,xnode,xgauss);
%         end
%         se = se + b*u;
%     end
    
    if mat.k2
        k2 = evalparam(mat,'k2',elem,xnode,xgauss);
        se = se + k2*(u.*u)*gradu;
    end
    
%     if mat.r
%         warning('mal calcule')
%         r = evalparam(mat,'r',elem,xnode,xgauss);
%         se = se + r/2*(u.*u);
%     end
    
%     if mat.r2
%         warning('mal calcule')
%         r2 = evalparam(mat,'r3',elem,xnode,xgauss);
%         se = se + r2/3*(u.*u.*u);
%     end
    
%     if mat.r3
%         warning('mal calcule')
%         r3 = evalparam(mat,'r3',elem,xnode,xgauss);
%         se = se + r3/4*(u.*u.*u.*u);
%     end
end


