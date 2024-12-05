function Z = integratemtimes(a,b,c)
% function Z = integratemtimes(a,b)
% calcul de integrale de (a*b)
%
% function Z = integratemtimes(a,b,c)
% calcul de integrale de (a*b*c)
% 
switch nargin
    case 2
     if ~istime(a)
      Z = a*integrate(b);    
     elseif ~istime(b)
      Z = integrate(a)*b;   
     elseif isa(a,'TIMEMATRIX') && isa(b,'TIMERADIALMATRIX')
       if all(size(a)==1)
        M = getMmatrix(a);     
        B = double(a)*M*double(b.L)';
        Z = reshape(double(b.V)*B',size(b));
       else
        M = getMmatrix(a);   
        B = MULTIMATRIX(double(a)*M*double(b.L)',size(a),[a.m,1]);
        Z = multisum(a.V*B);    
       end
     elseif isa(b,'TIMEMATRIX') && isa(a,'TIMERADIALMATRIX')
       if all(size(b)==1)
        M = getMmatrix(b);     
        B = double(b)*M*double(a.L)';
        Z = reshape(double(a.V)*B',size(a));
       else
        M = getMmatrix(b);   
        B = MULTIMATRIX(double(b)*M*double(a.L)',size(b),[a.m,1]);
        Z = multisum(a.V*B);    
       end
     elseif isa(b,'TIMERADIALMATRIX') && isa(a,'TIMERADIALMATRIX')
        if all(size(b)==1) 
        Z = integratemtimes(a,expand(b));            
        elseif all(size(a)==1)
        Z = integratemtimes(expand(a),b);   
        else
        LL = integratemtimes(a.L,b.L');
        B = MULTIMATRIX(double(b.V)*LL',size(b),[a.m,1]); 
        Z = multisum(a.V*B);
        end
     else
      error('pas defini')
     end
     
    case 3
        
        
    if ~istime(a)
      Z = a*integratemtimes(b,c);
    elseif ~istime(b)    
      error('pas programme')
    elseif ~istime(c)
      Z = integratemtimes(a,b)*c;    
    elseif isa(a,'TIMEMATRIX') && isa(b,'TIMEMATRIX') 
        error('a faire')
      
    elseif isa(b,'TIMEMATRIX') && isa(c,'TIMEMATRIX')
        error('a faire')
      
    elseif isa(a,'TIMEMATRIX') && isa(c,'TIMEMATRIX')
        error('a faire')
      
    elseif isa(a,'TIMEMATRIX')
      if all(size(a)==1)
         Z = integratemtimes(a*b,c); 
      else
         Z = integratemtimes(a,b*c); 
      end
    elseif isa(b,'TIMEMATRIX')
       Z = integratemtimes(a*b,c); 
    elseif isa(c,'TIMEMATRIX')
       Z = integratemtimes(a*b,c); 
        
    elseif isa(b,'TIMERADIALMATRIX') && isa(c,'TIMERADIALMATRIX') && isa(a,'TIMERADIALMATRIX')
       error('a faire')
       
    
    else
    error('pas defini')    
        
    end
  
    
    otherwise
      error('pas le bon nombre d''arguments')   

        
        
end

