function an = prodscal(a,b,M1,M2)

if isa(a,'cell') && isa(b,'cell')

    if length(a)~=length(b) || any(size(a{1})~=size(b{1}))
        error('ca correspond pas')
    end

    if nargin>=3 && ~isempty(M1) && size(M1,1)~=numem(a{1})
        error('la premiere metrique n''a pas la bonne dimension')
    end

    if nargin>=4 && ~isempty(M2) && length(a)~=size(M2,1)
        error('la deuxieme metrique n''a pas la bonne dimension')
    end


    switch nargin
    case 2
        an=0;       
        for k=1:length(a)
            an=an+prodscal(a{k},b{k});    
        end
    case 3
        an=0;       
        for k=1:length(a)
            an=an+prodscal(a{k},b{k},M1);    
        end        
    case 4
        if ~isempty(M1)  
            for i=1:length(a)
                a{i}=M1*a{i};    
            end  
            an = prodscal(a,b,[],M2);  

        elseif isempty(M2)
            an = prodscal(a,b);   

        else

            t1 = sqrt(sum(sum(M2.^2)));
            t2 = sqrt(sum(sum((M2-diag(diag(M2))).^2)));
            if t2/t1<1e-10
                for i=1:length(a)
                    a{i}=a{i}*M2(i,i);    
                end
            else
                c=cell(1,length(a));
                for k=1:length(c)
                    c{k}=a{1}*M2(k,1);
                    for l=2:length(a)
                        c{k} = c{k}+ a{l}*M2(k,l);    
                    end
                end
                a=c;
            end

            an = prodscal(a,b);   
        end



    otherwise 
        error('pas defini')
    end

else

    switch nargin
    case 3
        if ~isempty(M1) 
            if size(b,1)~=size(M1,2) 
                error('ca correspond pas')
            end
            b=M1*b;
        end
        an = prodscal(a,b);

    case 4
        if ~isempty(M2) 
            if size(a,2)~=size(M2,1)
                error('ca correspond pas')
            end
            a=a*M2;
        end
        an = prodscal(a,b,M1);

    case 2

        an = full(sum(sum(a.*b)));    %an = full(a(:)'*b(:));

    otherwise 
        error('pas defini')
    end
end
