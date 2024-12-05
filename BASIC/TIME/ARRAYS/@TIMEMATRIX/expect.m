function w = expect(a,b,c)
% function w = expect(a,b,c)

switch nargin
    case 1
        w=a;
        if isa(a.value,'cell')
            for i=1:length(a.value)
                w.value{i}=expect(a.value{i});
            end
        else
            w.value = expect(a.value);
        end
    case 2
        if ~isa(a,'TIMEMATRIX')
            w=b;
            if isa(b.value,'cell')
                for i=1:length(b.value)
                    w.value{i}=expect(a,b.value{i});
                end
            else
                error('pas programme')
            end
            
        elseif ~isa(b,'TIMEMATRIX')
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
