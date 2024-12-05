function W = myorth(W,j)
% function W = myorth(W,j)

    if size(W.F,1)>1
        error('must be rank one')
    end
    if nargin==1
        j = 1:size(W.F,2);
    end
    for k=1:length(j)
       W.F{j(k)} = mygram(W.F{j(k)});    
    end

    