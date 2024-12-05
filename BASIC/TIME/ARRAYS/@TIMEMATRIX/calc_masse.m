function Mt = calc_masse(a)

if numel(a)>1
    error('marche pour un scalaire')
end

dt = getdt(a);
nt = getnt(a);

switch getapproxparam(a,'type')
     case 'DG'
    
         switch  getapproxparam(a,'p')   
     
             case 0
         Mt = spdiags(dt(:).*(a.value(:)),0,nt,nt);        
             case 1
         
         %matelem = [2,1;1,2]/6;
         matelem{1} = [3,1;1,1]/12;
         matelem{2} = [1,1;1,3]/12;
         for i=1:2
             adt{i} = a.value(i:2:end);
             adt{i} = adt{i}(:).*dt(:);
         end
         
       
         I=[];
         J=[];
         V=[];
         for i=1:2
             for j=1:2
         I = [I;[i:2:nt*2]'];        
         J = [J;[j:2:nt*2]'];    
         Vtemp =adt{1}*matelem{1}(i,j)+adt{2}*matelem{2}(i,j);
         V = [V;Vtemp];
             end
         end
         Mt = sparse(I,J,V);

             case 2
         a = meanperinterval(a);
         dt = dt(:).*a(1:3:end);
         %matelem = [4,2,-1;2,16,2;-1,2,4]/30;
         matelem{1} = [ 39, 20, -3 ; ...
                      20, 16, -8;...
                      -3, -8, -3]/420;
         matelem{2} =  [ 20,  16,  -8;...
                      16, 192,  16;...
                      -8,  16,  20]/420;
         matelem{3} = [ -3, -8, -3;...
                      -8, 16, 20;...
                      -3, 20, 39]/420;      
         for i=1:3
             adt{i} = a.value(i:3:end);
             adt{i} = adt{i}(:).*dt(:);
         end
         
         I=[];
         J=[];
         V=[];
         for i=1:3
             for j=1:3
         I = [I;[i:3:nt*3]'];        
         J = [J;[j:3:nt*3]']; 
         Vtemp =adt{1}*matelem{1}(i,j)+adt{2}*matelem{2}(i,j)+...
               adt{3}*matelem{3}(i,j) ;
         V = [V;Vtemp];
             end
         end
         Mt = sparse(I,J,V);
             otherwise
         error('pas defini')

         end
     
    otherwise
        error('pas defini')
         
     
 end
        