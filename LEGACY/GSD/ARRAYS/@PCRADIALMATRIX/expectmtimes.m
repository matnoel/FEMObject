function Z = expectmtimes(a,b,c)
% function Z = expectmtimes(a,b)
% calcul de E(a*b)
%
% function Z = expectmtimes(a,b,c)
% calcul de E(a*b*c)

switch nargin
    case 2
     if ~israndom(a)
      Z = a*expect(b);    
     elseif ~israndom(b)
      Z = expect(a)*b;   
     elseif isa(a,'PCMATRIX') && isa(b,'PCRADIALMATRIX')
        
        
        if all(size(a)==1) && ~iscell(b)
        A = double(a)*double(b.L)';
        Z = reshape(sparsemtimes(double(b.V),A'),size(a));
        else
       
        a = MULTIMATRIX(double(a)*double(b.L)',size(a),[b.m,1]);
        Z = multisum(a*b.V);     
        end  
       
        
     elseif isa(b,'PCMATRIX') && isa(a,'PCRADIALMATRIX')


if all(size(b)==1) && ~iscell(a)
        B = double(b)*double(a.L)';
        Z = reshape(sparsemtimes(double(a.V),B'),size(a));
       else
        B = MULTIMATRIX(double(b)*double(a.L)',size(b),[a.m,1]);
        Z = multisum(a.V*B);    
       end
       
     elseif isa(b,'PCRADIALMATRIX') && isa(a,'PCRADIALMATRIX')
        if all(size(b)==1) 
        Z = expectmtimes(a,expand(b));            
        elseif all(size(a)==1)
        Z = expectmtimes(expand(a),b);   
        else
        LL = expectmtimes(a.L,b.L');
        if iscell(b.V)
        B = multimtimes(sparse(LL),b.V);    
        else
        B = MULTIMATRIX(double(b.V)*sparse(LL)',size(b),[a.m,1]); 
        end
        
        Z = multisum(a.V*B);
        end
     else
      error('pas defini')
     end
     
    case 3
      
    if ~israndom(a)
      Z = a*expectmtimes(b,c);
    elseif ~israndom(b)    
      Z = expectmtimes(a,mtimes(b,c)); 
    elseif ~israndom(c)
      Z = expectmtimes(a,b)*c;    
    elseif isa(a,'PCMATRIX') && isa(b,'PCMATRIX') 
      if all(size(a)==1) && all(size(b)==1)
         Z = expecttimes(a,b,c);
      elseif all(size(a)==1) %% a optimiser et traiter le cas size(b)==1
          Z = expectmtimes(b,a*c);
      else
          Z = expectmtimes(a,b*c);  
      end
    elseif isa(b,'PCMATRIX') && isa(c,'PCMATRIX')
      if all(size(b)==1) && all(size(c)==1)
          Z = expecttimes(a,b,c);
      elseif all(size(c)==1)
          Z = expectmtimes(a*c,b);
      else
          Z = expectmtimes(a*b,c);  
      end
    elseif isa(a,'PCMATRIX') && isa(c,'PCMATRIX')
      if all(size(a)==1) && all(size(c)==1)
          Z = expecttimes(a,b,c);
      elseif all(size(a)==1)
          Z = expectmtimes(a*b,c);
      else
          Z = expectmtimes(a,b*c);  
      end
    elseif isa(a,'PCMATRIX')
      if all(size(a)==1)
         Z = expectmtimes(a*b,c); 
      else
         Z = expectmtimes(a,b*c); 
      end
    elseif isa(b,'PCMATRIX')
      if all(size(b)==1)
         Z = expectmtimes(a,b*c); 
      else
         MAC = double(a.L)*getmasse(b)*double(c.L)';
         B = MULTIMATRIX(double(b)*double(MAC)',size(b),[a.m,c.m]);
         Z = multisum(a.V{1}*B{1,1:c.m}*c.V);
         for i=2:a.m
         Z = Z + multisum(a.V{i}*B{i,:},c.V);
         end
      end
    elseif isa(c,'PCMATRIX')
      if all(size(c)==1)
         Z = expectmtimes(a,b*c); 
      else
         Z = expectmtimes(a*b,c);   
      end
        
    elseif isa(b,'PCRADIALMATRIX') && isa(c,'PCRADIALMATRIX') && isa(a,'PCRADIALMATRIX')
        if iscalcmasse(a)
       
        MBC = double(b.L)*getmasse(a)*double(c.L)';
        if iscell(a.V)
        A = multimtimes(double(MBC),double(a.V));    
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
        A = multimtimes(MCB',a.V);    
        A = reshapem(A,[c.m,b.m]);
        else
        A = MULTIMATRIX(double(a.V)*MCB,size(a),[c.m,b.m]);
        end
        BC = multitranspose(b.V)*c.V;
        Z = multisum(A*BC);
        
        elseif iscalcmasse(c)
        MAB = double(a.L)*getmasse(c)*double(b.L)'; 
        if iscell(c.V)
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
        
        Z = expectmtimes(a,b,c);
        
        end
    
    else
    error('pas defini')    
        
    end
  
    
    otherwise
      error('pas le bon nombre d''arguments')   

        
        
end

