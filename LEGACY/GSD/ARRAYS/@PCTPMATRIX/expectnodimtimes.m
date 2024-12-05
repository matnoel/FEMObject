function x = expectnodimtimes(nodim,x,a,b)


if nargin==3 && ~israndom(x)
    x = x.*expectnodim(nodim,a);
elseif nargin==3 && ~israndom(a)
    x = expectnodim(nodim,x).*a;
elseif nargin==4 && ~israndom(x)
    x = x.*expectnoditimes(nodim,a,b);
    
elseif nargin==3 && isa(x,'PCTPMATRIX') && isa(a,'PCTPMATRIX')
    
    %% E_\nodim(ax)= x_nodim*a_nodim * prod_i\nodim E(a_i*x_i)
    
    %% experance E_i(a_i*x_i) sur les dimensions autres que nodim
     isranda = isranddim(a);
     israndx = isranddim(x);
    for i=setdiff(1:getnbdim(x),nodim)
     if isranda(i) && israndx(i) % E(x*a)
      x.phi{i} = a.phi{i}'*x.phi{i}; 
     elseif israndx(i) % E(x)*a
      Hi = mean(x.POLYCHAOSTP,i);   
      x.phi{i} = a.phi{i}*(Hi'*x.phi{i});    
     elseif isranda(i) % E(a)*x
      Hi = mean(x.POLYCHAOSTP,i);   
      x.phi{i} = x.phi{i}*(Hi'*a.phi{i});   
     else % a*x
      x.phi{i} = a.phi{i}*x.phi{i}; 
     end
    end

     %% multplication a_nodim*x_nodim selon la dimension nodim
     isranda = isranddim(a);
     israndx = isranddim(x);
     for i=nodim
     if isranda(i) && israndx(i) 
      if isximasse(a)   
      x.phi{i} = get_ximasse(a,i)*x.phi{i}; 
      elseif isximasse(x)   
      x.phi{i} = get_ximasse(x,i)*a.phi{i}; 
      else
      x.phi{i} = getmatrix(multimtimes(a.phi{i}',getmasseuni(x,i)))*x.phi{i};     
      end
     else
      x.phi{i} = x.phi{i}*a.phi{i};   
     end
     end
     
     
x.phi0 = x.phi0.*a.phi0;
     x.ximasse={};
    x.isranddim(setdiff(1:getnbdim(x),nodim))=0;
    
elseif nargin==4 && isa(x,'PCTPMATRIX') && isa(a,'PCTPMATRIX') && isa(b,'PCTPMATRIX')
    
   
    
    %% E_\nodim(xab)= x_nodim*a_nodim*b_nodim * prod_i\nodim E_i(x_i*a_i*b_i)
    
    %% experance E_i(x_i*a_i*b_i) sur les dimensions autres que nodim
    
    if ~isximasse(x) && ~isximasse(a) && ~isximasse(b)
        warning('pourrait etre optimise en calculant ximasse')
        x = calc_ximasse(x);
    end
    
     isranda = isranddim(a);
     israndx = isranddim(x);
     israndb = isranddim(b);

     for i=setdiff(1:getnbdim(x.POLYCHAOSTP),nodim)
     if isranda(i) && israndx(i) && israndb(i)
      if isximasse(x)   
      x.phi{i} = a.phi{i}'*get_ximasse(x,i)*b.phi{i}; 
      elseif isximasse(a)   
      x.phi{i} = x.phi{i}'*get_ximasse(a,i)*b.phi{i}; 
      else
      x.phi{i} = x.phi{i}'*get_ximasse(b,i)*a.phi{i};     
      end
     elseif israndx(i) && isranda(i)
      x.phi{i} = (a.phi{i}'*x.phi{i})*b.phi{i};    
     elseif israndx(i) && israndb(i)
      x.phi{i} = a.phi{i}'*(x.phi{i}'*b.phi{i});    
     elseif isranda(i) && israndb(i)
      x.phi{i} = x.phi{i}*(a.phi{i}'*b.phi{i});    
     else
       
      if isranda(i)
      x.phi{i} = (mean(x.POLYCHAOSTP,i)'*a.phi{i})*x.phi{i}*b.phi{i};   
      elseif israndb(i)
      x.phi{i} = a.phi{i}'*x.phi{i}*(mean(x.POLYCHAOSTP,i)'*b.phi{i});      
      elseif israndx(i)
      x.phi{i} = a.phi{i}'*(mean(x.POLYCHAOSTP,i)'*x.phi{i})*b.phi{i};     
      else
      x.phi{i} = a.phi{i}'*x.phi{i}*b.phi{i};    
      end
     end
    end
    
    % isranda = isranddim(a);
    % israndx = isranddim(x);
    % israndb = isranddim(b);
    for i=nodim
     if isranda(i) && israndx(i) && israndb(i)
      error('pas programme')  
     elseif israndx(i) && isranda(i)
      if isximasse(x)   
      x.phi{i} = (get_ximasse(x,i)*a.phi{i})*b.phi{i}; 
      else
      x.phi{i} = (get_ximasse(a,i)*x.phi{i})*b.phi{i};       
      end
     elseif israndx(i) && israndb(i)
      if isximasse(x)   
      x.phi{i} = (get_ximasse(x,i)*b.phi{i})*a.phi{i}; 
      else
      x.phi{i} = (get_ximasse(b,i)*x.phi{i})*a.phi{i};     
      end
     elseif isranda(i) && israndb(i)
      if isximasse(a)   
      x.phi{i} = (get_ximasse(a,i)*b.phi{i})*x.phi{i}; 
      else
      x.phi{i} = (get_ximasse(b,i)*a.phi{i})*x.phi{i};     
      end
     else
      x.phi{i} = x.phi{i}*a.phi{i}*b.phi{i};
     end
    end
    
    x.phi0 = x.phi0.*a.phi0.*b.phi0;
    x.ximasse={};
    x.isranddim(setdiff(1:getnbdim(x),nodim))=0;
    
else
    error('pas programme')
end

