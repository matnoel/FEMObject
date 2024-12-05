function Z = expect(A,b,c)
% function Z = expect(A,b)
% calcul de Z = E(A*b)
% avec Zi = E(A*b_i)
%
% function Z = expect(A,b,c)
% calcul de Z = E(A*b*c)
% avec Zij = E(A*b_i*c_j)

switch nargin
    case 1
        Z = mean(A);
    case 2
        if isempty(b)
            Z = expect(A);
            return
        end
        if ~israndom(A)
            if isa(A,'MYDOUBLEND')
                At = double(A);
                B = mean(b);
                Z = At(:) * B(:)';
                Z = reshape(full(Z),[size(A),prod(size(B))]);
                Z = MYDOUBLEND(Z);
                if prod(size(B))>1
                    sA = size(A);
                    Z = MULTIMYDOUBLEND(Z,length(sA)+1,size(B));
                end
            else
                B = sparse(mean(b));
                Z = A(:) * B(:)';
                Z = MULTIMATRIX(Z,size(A),size(B));
            end
        else
            if isa(b,'PCMATRIX')
                % Z = multimtimesold(A.MULTIMATRIX,double(b)');
                if iscell(b)
                    b = cell2mat(b);
                end
                if isa(A,'PCMATRIX')
                    Z = multimtimes(double(b),A.MULTIMATRIX);
                end
                
            else
                error('pas programme');
            end
        end
    case 3
        if isempty(c)
            Z = expect(A,b);
            return
        end
        
        if ~israndom(A)
            if isa(A,'MYDOUBLEND')
                At = double(A);
                B = sparse((double(b)*double(c)'));
                Z = At(:) * B(:)';
                Z = reshape(full(Z),[size(A),prod(size(B))]);
                Z = MYDOUBLEND(Z);
                if prod(size(B))>1
                    sA = size(A);
                    Z = MULTIMYDOUBLEND(Z,length(sA)+1,size(B));
                end
            else
                B = sparse((double(b)*double(c)'));
                Z = A(:) * B(:)';
                Z = MULTIMATRIX(Z,size(A),size(B));
            end
        else
            if isa(b,'PCMATRIX') && isa(c,'PCMATRIX')
                masse=getmasse(A);
                masseab = sparse(double(b))*masse*sparse(double(c)');
                
                % Z =  multimtimesold(A.MULTIMATRIX,double(masseab)');
                Z = multimtimes(double(masseab),A.MULTIMATRIX);
                Z  = reshapem(Z,[numel(b),numel(c)]);
            elseif isa(b,'PCMATRIX') && isempty(c)
                if isa(A,'PCMATRIX')
                    Z = multimtimes(double(b),A.MULTIMATRIX);
                end
            else
                
                
                error('pas programme');
            end
            
        end
        
end

if isa(Z,'MULTIMATRIX') && length(Z)==1
    Z=reshape(double(Z),size(A));
end
