function rep = polycmp(p1,p2)

if ~strcmp(p1.type,p2.type)
    rep=0;
else
    param1 = p1.param;
    param2 = p2.param;
    
    if strcmp(p1.type,'fe') || strcmp(p1.type,'wavelets')
        param1 = rmfield(param1,'p');
        param2 = rmfield(param2,'p');
    end
    
    
    %%%%%%%%%%%%ajout� par elias
    if strcmp(p1.type,'felagrange')
        param1 = rmfield(param1,'p');
        param2 = rmfield(param2,'p');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    t1=struct2cell(param1);
    t2=struct2cell(param2);
    rep=1;
    
    if length(t1)~=length(t2)
        rep=0;
    else
        for i=1:length(t1)
            %%%%%%%%%ajout� par Elias%%%%%%%%%%
            if ~strcmp(class(t1{i}),class(t2{i}))
                rep=0;break;
 
            elseif iscell(t1{i})&&iscell(t2{i})
                A=double(MYDOUBLEND(MULTIMATRIX(t1{i})));
                B=double(MYDOUBLEND(MULTIMATRIX(t2{i})));
                if ~(numel(t1{i})==numel(t2{i})) || ~all(all(~(A-B)))
                    rep=0;
                    break
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif isa(t1{i},'RANDVAR') && isa(t2{i},'RANDVAR')
                rep = (t1{i}==t2{i}); 
            else
                if ~(numel(t1{i})==numel(t2{i})) || ~all(all(~(t1{i}-t2{i})))
                    rep=0;
                    break
                end
                
            end
        end
    end
    
end