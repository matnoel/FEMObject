function [A,b,PC]=pcsystemupdate(A,b,varargin)
% function [A,b,PC]=pcsystemupdate(A,b)

if ischarin('noupdate',varargin)
    try
        PC = getPC(b);
    catch
        PC = getPC(A);
    end
else
    PC = getclassin('POLYCHAOS',varargin);
    if isempty(PC)
        if isa(b,'POLYCHAOS')
            PC = getPC(b);
        else
            PC = getPC(A);
        end
    else
        PC = getPC(PC);
        if isa(b,'POLYCHAOS') && ~polycmp(b,PC)
            if isin(b,PC)
                b=prolonge(b,PC);
            else
                b=project(b,PC);
            end
        end
    end
    
    if isa(b,'double')
        b=b*one(PC);
    end
    
    if israndom(A)
        A = calc_masse(A,PC);
    end
    
end