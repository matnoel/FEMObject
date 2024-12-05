function Z = integratetimes(a,b,c)
% function Z = integratetimes(a,b)
% calcul de integrale de (a.*b)
%
% function Z = integratetimes(a,b,c)
% calcul de integrale de (a.*b.*c)
% 
switch nargin
    case 2
     if ~istime(a)
      Z = a.*integrate(b);    
     elseif ~istime(b)
      Z = integrate(a).*b;   
     elseif isa(a,'TIMEMATRIX') && isa(b,'TIMERADIALMATRIX')
        Z = integratetimes(b,a);  
     elseif isa(b,'TIMEMATRIX') && isa(a,'TIMERADIALMATRIX')
       if all(size(b)==1)
        M = getMmatrix(b);     
        B = double(b)*M*double(a.L)';
        Z = reshape(double(a.V)*B',size(a));
       else
        M = getMmatrix(b);   
        B = MULTIMATRIX(double(b)*M*double(a.L)',size(b),[a.m,1]);
        Z = multisum(a.V.*B);    
       end
     elseif isa(b,'TIMERADIALMATRIX') && isa(a,'TIMERADIALMATRIX')
        if all(size(b)==1) 
        Z = integratetimes(a,expand(b));            
        elseif all(size(a)==1)
        Z = integratetimes(expand(a),b);   
        else
        LL = integratemtimes(a.L,b.L');
        B = MULTIMATRIX(double(b.V)*LL',size(b),[a.m,1]); 
        Z = multisum(a.V.*B);
        end
     else
      error('pas defini')
     end
     
    case 3
        
        
    if ~istime(a)
      Z = a.*integratetimes(b,c);
    elseif ~istime(b)    
      Z = integratetimes(a,c).*b;
    elseif ~istime(c)
      Z = integratetimes(a,b).*c;    
    elseif isa(a,'TIMEMATRIX') && isa(b,'TIMEMATRIX') 
        error('a faire')
      if all(size(a)==1) & all(size(b)==1)
         if iscalcmasse(c)
             MC = getmasse(c);
         else
             MC = getmasse(calc_masse(c))
         end
             MAB = double(a)*MC*double(b)';
             Z = reshape(double(c.V)*double(MAB)',size(c));
     
      elseif all(size(a)==1)
          Z = integratetimes(b,a.*c);
      else
          Z = integratetimes(a,b.*c);  
      end
    elseif isa(b,'TIMEMATRIX') && isa(c,'TIMEMATRIX')
        error('a faire')
      if all(size(b)==1) & all(size(c)==1)
          if iscalcmasse(a)
             MA = getmasse(a);
         else
             MA = getmasse(calc_masse(a));
          end
         
             MBC = double(b)*MA*double(c)';
             Z = reshape(double(a.V)*double(MBC)',size(a));
     
      elseif all(size(c)==1)
          Z = integratetimes(a.*c,b);
      else
          Z = integratetimes(a.*b,c);  
      end
    elseif isa(a,'TIMEMATRIX') && isa(c,'TIMEMATRIX')
        error('a faire')
      if all(size(a)==1) & all(size(c)==1)
          if iscalcmasse(b)
             MB = getmasse(b);
         else
             MB = getmasse(calc_masse(b));
         end
             MAC = double(a)*MB*double(c)';
            
             Z = reshape(double(b.V)*double(MAC)',size(b));
      elseif all(size(a)==1)
          Z = integratetimes(a.*b,c);
      else
          Z = integratetimes(a,b.*c);  
      end
    elseif isa(a,'TIMEMATRIX')
      if all(size(a)==1)
         Z = integratetimes(a.*b,c); 
      else
         Z = integratetimes(a,b.*c); 
      end
    elseif isa(b,'TIMEMATRIX')
       Z = integratetimes(b,a,c); 
    elseif isa(c,'TIMEMATRIX')
       Z = integratetimes(c,a,b); 
        
    elseif isa(b,'TIMERADIALMATRIX') && isa(c,'TIMERADIALMATRIX') && isa(a,'TIMERADIALMATRIX')
       error('a faire')
       if iscalcmasse(a)
       
        MBC = double(b.L)*getmasse(a)*double(c.L)';
        A = MULTIMATRIX(double(a.V)*double(MBC)',size(a),[b.m,c.m]);
        BC = b.V*multitranspose(c.V);
        Z = multisum(A*BC);
        
        elseif iscalcmasse(b)
        MAC = double(a.L)*getmasse(b)*double(c.L)';
        MCB = reshape(double(MAC),a.m,c.m*b.m);
        A = MULTIMATRIX(double(a.V)*MCB,size(a),[c.m,b.m]);
        BC = multitranspose(b.V)*c.V;
        Z = multisum(A*BC);
        
        elseif iscalcmasse(c)
        MAB = double(a.L)*getmasse(c)*double(b.L)'; 
        C = MULTIMATRIX(double(c.V)*double(MAB)',size(c),[a.m,b.m]);
        AB = a.V*multitranspose(b.V);
        Z = multisum(AB*C);    
        else
        if a.m <= b.m & a.m <=c.m 
          a=calc_masse(a);         
        elseif b.m <= a.m & b.m <=c.m 
        b=calc_masse(b);   
        elseif c.m <= a.m & c.m <=b.m 
        c=calc_masse(c);  
        end
        
        Z = integratetimes(a,b,c);
        
        end
    
    else
    error('pas defini')    
        
    end
  
    
    otherwise
      error('pas le bon nombre d''arguments')   

        
        
end

