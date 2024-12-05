function plot_lsdomainsol(elem,node,q,ls,varargin)
% function plot_lsdomainsol(elem,node,q,ls,varargin)

switch getlstype(elem)
    case {'in','indomain'}
        plot_sol(elem,node,q,varargin{:});
    case 'cut'
        connec = calc_conneclocal(elem,node);
        nodecoord = getcoord(node);
        lsval = getvalue(ls);
        xnode = getcoord(node,elem);
        
        issigma = ischarin('sigma',varargin);
        ampl = getcharin('ampl',varargin);
        varargin=delcharin('ampl',varargin);
        
        if issigma
            ksigma = getcharin('sigma',varargin);
            varargin = delcharin('sigma',varargin);
            varargin = setcharin('facecolor',varargin,'flat');
        end
        
        qelem=localize(elem,q);
        mat = getmaterial(elem);
        if issigma
            D = calc_opmat(mat,elem);
            B = calc_B(elem,xnode,[]);
        end
        
        for e=1:getnbelem(elem)
            % eleme=getelem(elem,e);
            connece = connec(e,:);
            nodecoorde = nodecoord(connece,:);
            
            if issigma
                Be = double(B(:,:,e));
            end
            
            lse = lsval(connece);
            
            qe = double(qelem(:,:,e));
            ae = reshape(qe,3,4);
            
            if all(lse<=0)
                
                if issigma
                    se = D*(Be*qe);
                    se = double(sigmacompo(se,ksigma,elem));
                    se = se(:);
                    varargin = setcharin('facevertexcdata',varargin,se);
                end
                
                patch('faces',1:3,'vertices',nodecoorde+ampl*ae',varargin{:});
                
            elseif all(lse>=0)
                
            else
                % division en sous-elements
                [connecin,connecout,xlnodeplus]=lsdivide_oneelem(lse);
                
                xlnodetotal=[nodelocalcoordtet4();xlnodeplus];
                % nodecoordeplus = calc_x(eleme,nodecoorde,xlnodeplus);
                nodecoordeplus = Ntet4(xlnodeplus)*nodecoorde;
                nodecoordeplus = [nodecoorde;nodecoordeplus];
                Ne = Ntet4(xlnodetotal);
                ue1 = Ne*ae';
                
                if issigma
                    se1 = D*(Be*qe);
                    se1 = double(sigmacompo(se1,ksigma,elem));
                    se1 = repmat(se1(:),size(connecin,1),1);
                    varargin = setcharin('facevertexcdata',varargin,se1);
                end
                
                if size(connecin,1)>0
                    patch('faces',connecin,'vertices',nodecoordeplus+ampl*ue1,varargin{:});
                end
            end
            
            
        end
        
end


function se = sigmacompo(se,ksigma,elem)

if isa(ksigma,'char') && strcmp(ksigma,'mises')
    switch getindim(elem)
        case 1
            se = double(se(1));
        case 2
            tracese = 1/3*(se(1)+se(2));
            se(1) = se(1) - tracese;
            se(2) = se(2) - tracese;
            se = sqrt(3/2*(se(1)*se(1) + se(2)*se(2) + 2*se(3)*se(3)));
        case 3
            tracese = 1/3*(se(1)+se(2)+se(3));
            se(1) = se(1) - tracese;
            se(2) = se(2) - tracese;
            se(3) = se(3) - tracese;
            se = sqrt(3/2*(se(1)*se(1) + se(2)*se(2) + se(3)*se(3)...
                + 2*(se(4)*se(4) + se(5)*se(5) + se(6)*se(6))));
    end
else
    se = se(ksigma);
end

return
