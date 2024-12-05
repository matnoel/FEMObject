function Z = expecttimes(a,b,c)
% function Z = expecttimes(a,b)
% calcul de E(a.*b)
%
% function Z = expecttimes(a,b,c)
% calcul de E(a.*b.*c)

switch nargin
    case 2
     if ~israndom(a)
      Z = a.*expect(b);    
     elseif ~israndom(b)
      Z = expect(a).*b;   
     elseif isa(a,'PCMATRIX') && isa(b,'PCRADIALMATRIX')
        Z = expecttimes(b,a);  
     elseif isa(b,'PCMATRIX') && isa(a,'PCRADIALMATRIX')
       if all(size(b)==1) && ~iscell(a)
        B = double(b)*double(a.L)';
        Z = reshape(double(a.V)*B',size(a));
       else
        B = MULTIMATRIX(double(b)*double(a.L)',size(b),[a.m,1]);
        Z = multisum(a.V.*B);    
       end
     elseif isa(b,'PCRADIALMATRIX') && isa(a,'PCRADIALMATRIX')
        if all(size(b)==1) 
        Z = expecttimes(a,expand(b));            
        elseif all(size(a)==1)
        Z = expecttimes(expand(a),b);   
        else
            
        LL = expectmtimes(a.L,b.L');
        if iscell(b.V)
        B = multimtimes(sparse(LL),b.V);
        else
        B = MULTIMATRIX(double(b.V)*LL',size(b),[a.m,1]); 
        end
        Z = multisum(a.V.*B);
        end
     else
      error('pas defini')
     end
     
    case 3
        
        
    if ~israndom(a)
      Z = a.*expecttimes(b,c);
    elseif ~israndom(b)    
      Z = expecttimes(a,c).*b;
    elseif ~israndom(c)
      Z = expecttimes(a,b).*c;    
    elseif isa(a,'PCMATRIX') && isa(b,'PCMATRIX') 
      if all(size(a)==1) && all(size(b)==1)
         if iscalcmasse(c)
             MC = getmasse(c);
         else
             MC = getmasse(calc_masse(c));
         end
         
         
             MAB = double(a)*MC*double(b)';
         
             if iscell(c.V)
             Z = multimtimes(double(MAB),c.V);
             Z = Z{1};
             else
             Z = reshape(double(c.V)*double(MAB)',size(c));
             end
             
      elseif all(size(a)==1)
          Z = expecttimes(b,a.*c);
      else
          Z = expecttimes(a,b.*c);  
      end
    elseif isa(b,'PCMATRIX') && isa(c,'PCMATRIX')
      if all(size(b)==1) && all(size(c)==1)
          if iscalcmasse(a)
             MA = getmasse(a);
         else
             MA = getmasse(calc_masse(a));
          end
         
             MBC = double(b)*MA*double(c)';
             
             if iscell(a.V)
             Z = multimtimes(double(MBC),a.V);
             Z = Z{1};
             else
             Z = reshape(double(a.V)*double(MBC)',size(a));
             end
             
      elseif all(size(c)==1)
          Z = expecttimes(a.*c,b);
      else
          Z = expecttimes(a.*b,c);  
      end

    elseif isa(a,'PCMATRIX') && isa(c,'PCMATRIX')
      if all(size(a)==1) && all(size(c)==1)
          if iscalcmasse(b)
             MB = getmasse(b);
         else
             MB = getmasse(calc_masse(b));
         end
             MAC = double(a)*MB*double(c)';
             if iscell(b.V)
             Z = multimtimes(double(MAC),b.V);
             Z = Z{1};
             else
             Z = reshape(double(b.V)*double(MAC)',size(b));
             end
      elseif all(size(a)==1)
          Z = expecttimes(a.*b,c);
      else
          Z = expecttimes(a,b.*c);  
      end
    elseif isa(a,'PCMATRIX')
      if all(size(a)==1)
         Z = expecttimes(a.*b,c); 
      else
         Z = expecttimes(a,b.*c); 
      end
    elseif isa(b,'PCMATRIX')
       Z = expecttimes(b,a,c); 
    elseif isa(c,'PCMATRIX')
       Z = expecttimes(c,a,b); 
        
    elseif isa(b,'PCRADIALMATRIX') && isa(c,'PCRADIALMATRIX') && isa(a,'PCRADIALMATRIX')
        if iscalcmasse(a)
       
        MBC = double(b.L)*getmasse(a)*double(c.L)';
        if iscell(a.V)
        A = multimtimes(double(MBC),a.V);
        A = reshapem(A,[b.m,c.m]);
        else
        A = MULTIMATRIX(double(a.V)*double(MBC)',size(a),[b.m,c.m]);
        end
        BC = b.V*multitranspose(c.V);
        Z = multisum(A*BC);
        
        elseif iscalcmasse(b)
        MAC = double(a.L)*getmasse(b)*double(c.L)';
        MCB = reshape(double(MAC),a.m,c.m*b.m);
        if iscell(a.V)
        A = multimtimes(double(MCB)',a.V);
        A = reshapem(A,[c.m,b.m]);
        else
        A = MULTIMATRIX(double(a.V)*MCB,size(a),[c.m,b.m]);
        end
        BC = multitranspose(b.V)*c.V;
        Z = multisum(A*BC);
        
        elseif iscalcmasse(c)
        MAB = double(a.L)*getmasse(c)*double(b.L)'; 
        if iscell(a.V)
        C = multimtimes(double(MAB),c.V);
        C = reshapem(C,[a.m,b.m]);
        else
        C = MULTIMATRIX(double(c.V)*double(MAB)',size(c),[a.m,b.m]);
        end
        AB = a.V*multitranspose(b.V);
        Z = multisum(AB*C);    
        else
        if a.m <= b.m && a.m <=c.m 
          a=calc_masse(a);         
        elseif b.m <= a.m && b.m <=c.m 
        b=calc_masse(b);   
        elseif c.m <= a.m && c.m <=b.m 
        c=calc_masse(c);  
        end
        
        Z = expecttimes(a,b,c);
        
        end
    
    else
    error('pas defini')    
        
    end
  
    
    otherwise
      error('pas le bon nombre d''arguments')   

        
        
end

