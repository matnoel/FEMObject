function a = convertradial(a)
% function u = convertradial(a)
% a : PCMATRIX
% u : PCRADIALMATRIX


    if numel(a)==1
    a = PCRADIALMATRIX(1,[1,1],a);
    else
    a = reshape(a,[numel(a),1]);
    V = diag(ones(1,numel(a)));
    a = PCRADIALMATRIX(V,size(a),a);    
    try 
    a = setmasse(a,getximasse(a));
    end
        
    end
 
    