function D = calc_opmatpc(mat,elem,xnode,xgauss,PC)
% function D = calc_opmatpc(mat,elem,xnode,xgauss,PC)

if nargin<=5
    PC = [];
end

switch getdim(elem)
    case 1
        E = evalparampc(mat,'E',PC,elem,xnode,xgauss);
        S = evalparampc(mat,'S',PC,elem,xnode,xgauss);
        D = E*S;
    case 2
        if isaxi(elem)
            E = evalparampc(mat,'E',PC,elem,xnode,xgauss);
            nu = evalparampc(mat,'NU',PC,elem,xnode,xgauss);
            D = E/(1+nu)/(1-2*nu)*[(1-nu),nu,nu,0;nu,(1-nu),nu,0;nu,nu,(1-nu),0;0,0,0,(1-2*nu)/2];
        else
            if israndom(getparam(mat,'DIM3')) || israndom(getparam(mat,'nu'))
                error('pas programme')
            end
            E = evalparampc(mat,'E',PC,elem,xnode,xgauss);
            e = evalparampc(mat,'DIM3',PC,elem,xnode,xgauss);
            nu = evalparampc(mat,'NU',PC,elem,xnode,xgauss);
            switch getoption(elem)
                case 'DEFO'
                    D = e*E/(1+nu)/(1-2*nu)*[(1-nu),nu,0;nu,(1-nu),0;0,0,(1-2*nu)/2];
                otherwise
                    
                    D = e*E/(1-nu^2)*[1,nu,0;nu,1,0;0,0,(1-nu)/2];
                    
            end
        end
    case 3
        E = evalparampc(mat,'E',PC,elem,xnode,xgauss);
        nu = evalparampc(mat,'NU',PC,elem,xnode,xgauss);
        
        D = E/(1+nu)/(1-2*nu)*[(1-nu),nu,nu,0,0,0;...
            nu,(1-nu),nu,0,0,0;...
            nu,nu,(1-nu),0,0,0;...
            0,0,0,(1-2*nu)/2,0,0;...
            0,0,0,0,(1-2*nu)/2,0;...
            0,0,0,0,0,(1-2*nu)/2];
        
end


