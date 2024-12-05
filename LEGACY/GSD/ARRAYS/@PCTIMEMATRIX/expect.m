function w = expect(a,b,c)
% function w = expect(a,b,c)

switch nargin
    case 1
        if isa(a.value,'cell')
            for i=1:length(a.value)
                value{i}=expect(a.value{i});
            end
        else
            value  = expect(a.value);
        end
        
        w = TIMEMATRIX(value,a.TIMEMODEL,a.s);
    case 2
        if ~isa(a,'PCTIMEMATRIX')
            w=b;
            if isa(b.value,'cell')
                for i=1:length(b.value)
                    w.value{i}=expect(a,b.value{i});
                end
            else
                error('pas programme')
            end
            
        elseif ~isa(b,'PCTIMEMATRIX')
            w=a;
            if isa(a.value,'cell')
                for i=1:length(a.value)
                    w.value{i}=expect(a.value{i},b);
                end
                
            else
                error('pas programme')
            end
        end
        
        
    otherwise
        error('pas programme')
end
