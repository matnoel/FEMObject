function Mt = getMmatrix(L)
% function Mt = getMmatrix(L)

dt = getdt(L);
nt = getnt(L);

switch getapproxparam(L,'type')
    case {'default','CG'}
        switch  getapproxparam(L,'p')
            case 0
                Mt = spdiags(dt(:),0,nt,nt);
            case 1
                dtd=[dt/6,0]';
                dtu=[0,dt/6]';
                dtm = ([dt/3,0]+[0,dt/3])';
                Mt = spdiags([dtd,dtm,dtu],[-1,0,1],nt+1,nt+1);
                %dtm = ([dt/2,0]+[0,dt/2])';
                %Mt = spdiags(dtm,[0],nt+1,nt+1);
            otherwise
                error('pas defini')
        end
    case 'DG'
        switch  getapproxparam(L,'p')
            case 0
                Mt = spdiags(dt(:),0,nt,nt);
            case 1
                matelem = [2,1;1,2]/6;
                I=[];
                J=[];
                V=[];
                for i=1:2
                    for j=1:2
                        I = [I;[i:2:nt*2]'];
                        J = [J;[j:2:nt*2]'];
                        V = [V;dt'*matelem(i,j)];
                    end
                end
                Mt = sparse(I,J,V);
            case 2
                matelem = [4,2,-1;2,16,2;-1,2,4]/30;
                I=[];
                J=[];
                V=[];
                for i=1:3
                    for j=1:3
                        I = [I;[i:3:nt*3]'];
                        J = [J;[j:3:nt*3]'];
                        V = [V;dt'*matelem(i,j)];
                    end
                end
                Mt = sparse(I,J,V);
            otherwise
                error('pas defini')
        end
    otherwise
        error('pas defini')
end
