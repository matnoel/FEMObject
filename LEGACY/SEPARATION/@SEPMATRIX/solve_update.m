function [u,utilde] = solve_update(A,b,u,utilde,Ametric,param)
% function [u,utilde] = solve_update(A,b,u,utilde,Ametric,param)

n=size(b);
if param.updatedimtest
    param.updatedim = param.updatedim(n(param.updatedim)>=getm(u));
end

if ~isempty(param.updatedim)
    for kkk=1:param.update
        if param.display;fprintf(' update #%d, dim ',kkk);end
        u0 = u;
        W = gathervectors(u);
        W.alpha = 1;
        if param.adjoint==0 || ~param.updateadjoint
            Wtilde = W;
        else
            Wtilde = gathervectors(utilde);
            Wtilde.alpha = 1;
        end
        
        if param.ortho
            W = myorth(W,param.orthodim);
            Wtilde = myorth(Wtilde,param.orthodim);
        end
        
        WAW = Wtilde'*A*W;
        Wb = Wtilde'*b;
        
        for jj=1:length(param.updatedim)
            j = param.updatedim(jj);
            if param.display;fprintf('#%d ',j);end
            WAW.F(:,j)= A.F(:,j);
            Wb.F(:,j)=b.F(:,j);
            
            if ismember(j,param.updatedimapprox)
                M = reduceto2Dtensor(WAW,j);
                S = reduceto2Dtensor(Wb,j);
                solver = SEPSOLVER(getdim(M), 'tol',1e-5,'maxorder',max(size(M)),...
                    'update',1,'updatedim',2,'display',0);
                % nn = size(W.F{j},2);
                % Wj0 = splitvectors(SEPMATRIX({W.F{j},speye(nn)},ones(1,nn)));
                % Wj = solve(M,S-M*Wj0,solver);
                % W.F{j} = expand(Wj0 + Wj);
                Wj = solve(M,S,solver);
                W.F{j} = expand(Wj);
            else
                M = assembleblock(WAW,j);
                if param.updateeps>0
                    M = M+speye(size(M,1))*param.updateeps;
                end
                S = assembleblock(Wb,j);
                Wj = reshape(M\S,n(j),u.m);
                W.F{1,j}=Wj;
            end
            
            %%%
            % if j==1
            %     result.W{u.m}(kkk) = subspaceangle(double(getV(param.reference)),Wj);
            % end
            %%%
            
            if param.adjoint && param.updateadjoint
                u = splitvectors(W);
                u.alpha = ones(1,u.m);
                u = normalizefuns(u);
                AT=A';
                btemp = Ametric*u;
                Wb = W'*btemp;
                WAW = W'*AT*Wtilde;
                WAW.F(:,j)= AT.F(:,j);
                Wb.F(:,j)=btemp.F(:,j);
                M = assembleblock(WAW,j);
                M = M+speye(size(M,1))*param.updateeps;
                S = assembleblock(Wb,j);
                Wj = reshape(M\S,n(j),u.m);
                Wtilde.F{1,j}=Wj;
            else
                Wtilde = W;
            end
            
            WAW.F(:,j)= cellfun(@(C) Wtilde.F{j}'*C*W.F{j},A.F(:,j),'UniformOutput',0);
            Wb.F(:,j)= cellfun(@(C) Wtilde.F{j}'*C,b.F(:,j),'UniformOutput',0);
            
        end
        
        u = splitvectors(W);
        u.alpha = ones(1,u.m);
        u = normalizefuns(u);
        if param.adjoint
            utilde = splitvectors(Wtilde);
            utilde.alpha = ones(1,utilde.m);
            utilde = normalizefuns(utilde);
        end
        
        if kkk<param.update
            errup0= sqrt(min(u.alpha)^2/sum(u.alpha.^2));
            stagn = full(norm(u-u0)/norm(u0));
            if param.display;fprintf('- stagnation = %d\n',stagn); end
            if stagn<errup0*param.itercritupdate
                break
            end
        else
            if param.display;fprintf('\n'); end
        end
        
        
    end
end


function R = reduceto2Dtensor(A,j)

Z = A.F(:,j);
A.F(:,j) = [];
% n = size(Z{1},1);
% nb = size(A.F{1},1);
B = cell(A.m,1);
B(:) = {1};
for k=1:size(A.F,2)
    B = cellfun(@times,B,A.F(:,k),'UniformOutput',false);
end
R = SEPMATRIX([Z,B],A.alpha);

return

