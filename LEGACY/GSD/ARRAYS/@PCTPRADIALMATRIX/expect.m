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
        if ~(isa(b,'PCTPMATRIX') || isa(b,'PCTPMATRIXSUM'))
            error('pas programme')
        end
        if ~israndom(A)
            Z = A*expect(b);
            return
        elseif ~israndom(b)
            Z = expect(A)*b;
        end
        Lb = A.L;
        for i=1:length(Lb)
            Lb{i}= expect(Lb{i},b);
        end
        
        if numel(b)==1
            Z = multimtimes([Lb{:}],A.V);
        end
        
    case 3
        if isempty(c)
            Z = expect(A,b);
            return
        end
        
        if ~(isa(b,'PCTPMATRIX') || isa(b,'PCTPMATRIXSUM')) || ~(isa(c,'PCTPMATRIX') || isa(c,'PCTPMATRIXSUM'))
            error('pas programme')
        end
        
        if ~israndom(A)
            Z = A*expect(b,c);
            return
        end
        
        Lbc = A.L ;
        
        for i=1:length(Lbc)
            Lbc{i}= expect(Lbc{i},b,c);
        end
        
        if numel(b)==1 && numel(c)==1
            Z = multimtimes([Lbc{:}],A.V);
        end
        
end

if isa(Z,'MULTIMATRIX') && numelm(Z)==1
    if iscell(Z)
        Z = reshape(Z{1},size(A));
    else
        Z=reshape(double(Z),size(A));
    end
end


