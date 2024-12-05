function Z = expect(A,b,c)
% function Z = expect(A,b,c)
% calcul de E(A*b)
% calcul de E(A*b*c)


switch nargin
    case 1
        Z = mean(A);
    case 2
        if isempty(b)
            Z = expect(A);
            return
        end
        if ~isa(b,'PCMATRIX')
            if isa(b,'PCRADIALMATRIX')
                b = expand(b);
            else
                error('pas programme')
            end
        end
        if ~israndom(A)
            Z = expect(A,b);
            return
        end
        Lb= sparse(expect(A.L,b));
        % Z = multimtimesold(A.V , A.D * double(Lb));
        
        Lb = (A.D *double(Lb));
        if size(Lb,2)>1 || iscell(A.V)
            Z = multimtimes(Lb',A.V);
            Z = reshapem(Z,size(b));
        else
            Z = reshape(double(A.V)*Lb,size(A.V));
        end
        
    case 3
        if isempty(c)
            Z = expect(A,b);
            return
        end
        
        if ~isa(b,'PCMATRIX') || ~isa(c,'PCMATRIX')
            if isa(b,'PCRADIALMATRIX')
                b = expand(b);
            else
                error('pas programme')
            end
            if isa(c,'PCRADIALMATRIX')
                c = expand(c);
            else
                error('pas programme')
            end
        end
        if ~israndom(A)
            Z = expect(A,b,c);
            return
        end
        
        A = actualise_masse(A);
        Lbc = sparse(double(b) * (A.DLmasse * double(c)')) ;
        
        
        % Z = multimtimesold(A.V , double(Lbc)');
        
        if all(size(Lbc)>1) || iscell(A.V)
            Z = multimtimes(double(Lbc),A.V);
            Z = reshapem(Z,numel(b),numel(c));
        else
            Z = reshape(double(A.V)*double(Lbc)',size(A.V));
        end
        
end

if isa(Z,'MULTIMATRIX') && numelm(Z)==1
    if iscell(Z)
        Z = reshape(Z{1},size(A));
    else
        Z=reshape(double(Z),size(A));
    end
end


