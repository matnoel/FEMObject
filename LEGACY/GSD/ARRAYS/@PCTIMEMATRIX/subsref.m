function u = subsref(u,s)
% function u = subsref(u,s)

if length(s)==2
    u = subsref(u,s(1));
    u = subsref(u,s(2));
    
else
    switch s.type
        case '{}'
            if isa(u.value,'cell')
                u = subsref(u.value,s);
            elseif isa(u.value,'double') || israndom(u.value)
                u = reshape(u.value(:,s.subs{1}),u.s);
            elseif isa(u.value,'FEELEMFIELD')
                u = u.value(:,s.subs{1});
            else
                error('pas programme')
            end
            
            
        case '()'
            
            
            if isa(u.value,'double')
                warning('le champ devrait etre aleatoire')
                v = MULTIMATRIX(u.value,u.s,[length(u.TIMEMODEL),1]);
                v = subsref(v,s);
                u.s = size(v);
                u.value = double(v);
            elseif isa(u.value,'MULTIMATRIX') || isa(u.value,'PCMATRIX')
                if strcmp(u.value,'MULTIMATRIX')
                    warning('le champ devrait etre aleatoire')
                end
                m = length(u.value);
                U = MULTIMATRIX(double(u.value),u.s,[length(u.TIMEMODEL),m]);
                U = subsref(U,s);
                u.s = size(U);
                U = MULTIMATRIX(double(U),[prod(u.s),length(u.TIMEMODEL)],[m,1]);
                if isa(u.value,'PCMATRIX')
                    U = PCMATRIX(U,size(U),getPC(u.value));
                end
                u.value = U;
                
            elseif isa(u.value,'PCRADIALMATRIX')
                U = getV(u.value);
                m = length(U);
                U = MULTIMATRIX(double(U),u.s,[length(u.TIMEMODEL),m]);
                U = subsref(U,s);
                u.s = size(U);
                U = MULTIMATRIX(double(U),[prod(u.s),length(u.TIMEMODEL)],[m,1]);
                u.value = setV(u.value,U);
                
            elseif isa(u.value,'FEELEMFIELD')
                
                if length(s.subs)==1
                    s.subs=[s.subs {':'}];
                    u.value = subsref(u.value,s);
                    u.s(1) = size(u.value,1);
                else
                    warning('ambiguite sur le subsref')
                    u = subsref(u.value,s);
                end
                
                
            elseif  israndom(u.value)
                
                if (length(s.subs)==1 || isa(s.subs{2},'char'))  && u.s(2)==1
                    s.subs{2}=':';
                    if isa(u.value,'cell')
                        for k =1:length(u.value)
                            u.value{k} = subsref(u.value{k},s);
                        end
                        u.s(1)=size(u.value{1},1);
                    else
                        u.value = subsref(u.value,s);
                        u.s(1)=size(u.value,1);
                    end
                else
                    error('pas programme')
                end
            else
                error('pas normal')
            end
    end
    
end