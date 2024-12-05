function s = simplify(s)

if isempty(getdim(s.PRODUCTSPACE))
    s = s.factor;
else
    repsimp = zeros(1,length(s.phi));
    for i=1:length(s.phi)
    if isa(s.phi{i},'double') && numel(s.phi{i})==1
        repsimp(i)=1;
        s.factor = s.factor * s.phi{i};
    end
    end
    if all(repsimp)
    s = s.factor;
    elseif any(repsimp)
    remdim = getdim(s.PRODUCTSPACE,find(repsimp));    
    s.PRODUCTSPACE=removedim(s.PRODUCTSPACE,remdim);
    s.phi(find(repsimp))=[];
    end
        
end
